% Script written 1.10.2021 by Ewa Beldzik 
% Follow-up analyses - omissions

clear; clc

load('E:\SEF2\beh\beh_sef.mat');
load('E:\SEF2\eeg\matfiles\ic_mf_stim.mat')

%%% EEG time-frequncy analysis
% dataset parameters
sb_in = setdiff(1:37,[5,25]);
N = length(sb_in);
ntrials = 1120;
% behavioural data
T = T(ismember(T.Subject,sb_in),:);

% baseline period
% bIdx = EEG.times > -1000 & EEG.times <-900; 

% stimulus-locked time epoch
timeEEG = -1200:20:1200;
timeStimIdx = dsearchn(EEG.times',timeEEG'); 

%%% wavelet parameters
num_frex = 60;
range_frex = [ 2 30 ];
range_cycl = [ 2 7 ];
%logarithmically spaced frequency
frex  = logspace(log10(range_frex(1)),log10(range_frex(2)),num_frex);
nCycs = logspace(log10(range_cycl(1)),log10(range_cycl(2)),num_frex); 

% initialize output time-frequency data
tf = [];

for sb=sb_in
    
    dataY = ic_tc(sb).data;
    
    % FFT parameters
    wave_time = -2:1/EEG.srate:2;
    half_wave = (length(wave_time)-1)/2;
    nWave = length(wave_time);
    nData = size(dataY,1)*size(dataY,2);
    nConv = nWave+nData-1;     
        
    dataX = fft( reshape(dataY,1,[]),nConv);
    tf_temp = zeros(num_frex,size(dataY,1),size(dataY,2));
    bf_temp = zeros(num_frex,size(dataY,2));
    
    % loop over frequencies
    for fi = 1:num_frex

        % create wavelet and get its FFT
        s = nCycs(fi)/(2*pi*frex(fi));
        wavelet  = exp(2*1i*pi*frex(fi).*wave_time) .* exp(-wave_time.^2./(2*s^2));
        waveletX = fft(wavelet,nConv);
        waveletX = waveletX./max(waveletX);

        % run convolution
        as = ifft(waveletX.*dataX,nConv);
        as = as(half_wave+1:end-half_wave);
        as = reshape(as,size(dataY,1),size(dataY,2));
                              
        % compute power
        tf_temp(fi,:,:) = abs(as.^2);
        
    end
    % stimulus-locked data in first column, response-locked data in second
    tf(sb).data = tf_temp(:,timeStimIdx,:);    

end

%% average and plot TF maps
figure(1); clf;clc

types = {'hit','omi'};
thetaIdx = frex > 4 & frex < 13;
timeIdx = timeEEG>-1000 & timeEEG<0;

% prepare output
tfTyp = []; tfAll = []; 
thetaAve =[]; thetaMin = []; 

% averaged TR maps
for sb=sb_in  
    Tsb=T(T.Subject==sb,:); 

    %for each trials type           
    for trialtyp=1:length(types)  
        typIdx = (Tsb.type=='con' | Tsb.type=='incon')  & Tsb.accur==types(trialtyp)  & Tsb.badep == 1;
        tfTyp(trialtyp).data(sb,:,:) = 10*log10(mean(tf(sb).data(:,:,typIdx),3)); %

        % extract mean and min theta values
        thetaAve(sb,trialtyp) = mean(mean(tfTyp(trialtyp).data(sb,thetaIdx,timeIdx),2));
        how_many(sb,trialtyp) = sum(typIdx);
    end  
  
end

subplot(1,2,1)
diffmap = squeeze(mean(tfTyp(2).data(sb_in,:,:))) - squeeze(mean(tfTyp(1).data(sb_in,:,:)));
contourf(timeEEG/1000,frex,diffmap,50,'linecolor','none'); 
hold on
rectangle('Position', [-1 4 1 8],'LineStyle','--')

xlabel('time [ms]'); ylabel('frequency [Hz]'); title('laps - correct');
set(gca,'clim',[-1.2 1.2],'ydir','normal'); xticks(-1:0.5:1.2); colormap(jet(256)); % ,
set(gca,'fontsize', 11)

%
subplot(1,4,3); 
wdth = 0.3; % width of boxplot
colors = {[.3 .4 1] [1 .3 .2]};

for i = 1:length(types)
    X = thetaAve(sb_in,i);
    
    % jitter for raindrops
    jit = (rand(size(X)) - 0.5) * wdth + i;
    scatter(jit,X,'SizeData',10,'MarkerFaceColor',colors{i},'MarkerEdgeColor','none');
    hold on

end
hold on
boxplot(thetaAve(sb_in,:),'Widths',1.5*wdth, 'Colors',[0 0 0])

xlim([0.3 2*length(types)-1.3]); xticks(1:length(types));
xticklabels({'cor','laps'}); set(gca,'fontsize', 11)
ylim([-19 -2]);

[h, p] = ttest(thetaAve(sb_in,1),thetaAve(sb_in,2))


%% fMRI analysis
% load fMRI data
mypath  = 'E:\SEF2\mri\ica\res_frontal25\'; 
load(strcat(mypath,'labels.mat'));
ic_list=[15,25,18,14,21,16,5,22];

% list of task sessions
ses_list={'simon20','simon50','stroop20','stroop50'};

%%% set parameters
% define epoch for hemodynamic response
t = -5:0.1:9.9;
% define matrix for hdr & its peak
hdr = nan(height(T),length(t),length(ic_list));

for sbs=1:35    
    sb=sb_in(sbs);
    
    for ses=1:length(ses_list)
        % Loading ic timecourse
        load(strcat(mypath,'sef_itc_sb',num2str(sbs),'_ses',num2str(ses),'.mat'));
        itc(end+1:end+120,:)=NaN;

        SesIdx = T.Subject==sb & T.session==ses_list(ses); 
        Tsb = T(SesIdx,:); 
        ampses = nan(height(Tsb),length(ic_list));
        
        % make list of trials in each session
        tr_list=round((Tsb.stimOnsetTime-Tsb.firstfixOnsetTime(1))/100,0);           
        hdrSes=nan(length(tr_list),length(t),length(ic_list));
        
        %extract HDR of each trial
        for tr=1:length(tr_list)  
            if tr_list(tr)-sum(t<0)+1 > 0
            hdrSes(tr,:,:)=itc(tr_list(tr)-sum(t<0)+1:tr_list(tr)+sum(t>0)+1,ic_list);% - itc(tr_list(tr),ic_list);
            else
            hdrSes(tr,:,:)= NaN;
            end
        end
        
        hdr(SesIdx,:,:) = hdrSes;
         
    end

end

%% average and plot fMRI HDR
hFig = figure(2); clf;
set( hFig, 'units','normalized','outerposition',[.3 .3 .3 .6]);

hdrAve=[]; hdrTyp = []; errBar=[];
colors = {'b','r'};
how_many =[];

% average HDR for each subject and trial type
for trialtyp=1:length(types) 
    for sb=sb_in    
        typIdx = T.Subject==sb & (T.type=='con' | T.type=='incon') & T.accur==types(trialtyp) & T.badep==1; % 
        how_many(sb,trialtyp) = sum(typIdx);
        hdrTyp(trialtyp).data(sb,:,:) = nanmean(hdr(typIdx,:,:));
    end
end
pval = nan(length(t),1);

ic = find(ismember(labels(ic_list),'area9')); 
    

    [~,pval] = ttest(hdrTyp(1).data(:,:,ic),hdrTyp(2).data(:,:,ic));  
    pcor = 100*fdr_bh(pval,0.05); pcor(pcor==0)=NaN;  
    
    area(t,pcor,-100,'LineStyle','none','FaceColor',[.87 .87 .87]);
    hold on
    for trialtyp=1:length(types)
        hdrAve(trialtyp,:) = nanmean(hdrTyp(trialtyp).data(sb_in,:,ic));
        errBar(trialtyp,:) = nanstd(hdrTyp(trialtyp).data(sb_in,:,ic))/sqrt(length(sb_in));     
        shadedErrorBar(t,hdrAve(trialtyp,:),errBar(trialtyp,:),'lineprops',colors{trialtyp}); 
    end
    
    ylim([min(min(hdrAve))-0.2 max(max(hdrAve))+0.2]); xlim ([t(1) t(end)]); %*4
    title(strcat('omi -',labels(ic_list(ic))))

xlabel('time [s]');
set(gca,'fontsize', 13)

