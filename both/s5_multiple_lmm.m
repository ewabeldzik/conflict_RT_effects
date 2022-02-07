% Script written 1.02.2021 by Ewa Beldzik 
% coupling EEG & fMRI usuing 'multiple' LME


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
bIdx = EEG.times > -500 & EEG.times <-200; 

% stimulus-locked time epoch
time(1).t = 0:25:1000;
timeStimIdx = dsearchn(EEG.times',time(1).t'); 

% response-locked time epoch
time(2).t = -400:25:300;

%%% wavelet parameters
num_frex = 30;
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
        tf_temp(fi,:,:) = abs(as.^2)- mean(abs(as(bIdx,:).^2));

    end
    % stimulus-locked data in first column, response-locked data in second
    tf(sb,1).data = tf_temp(:,timeStimIdx,:);
    
    
    % prepare resp-locked tf maps
    Tsb = T(T.Subject==sb,:);  
    typIdx = Tsb.stimRT~=0;     
    tf(sb,2).data = zeros(num_frex,length(time(2).t),ntrials);   
    
    for tr=find(typIdx)'
        time2resp = time(2).t + Tsb.stimRT(tr);
        time2respIdx = dsearchn(EEG.times',time2resp');
        %tf(sb,2).data(:,tp,tr) = mean(tf_temp(:,time2respIdx,tr);  
        
        % downsampling ft maps by averaging 50 ms time windows
        for tp = 1:length(time(2).t)
            tf(sb,2).data(:,tp,tr) = mean(tf_temp(:,time2respIdx(tp)-6:time2respIdx(tp)+6,tr),2);
        end
    end

end

%% fMRI analysis
% load fMRI data
mypath  = 'E:\SEF2\mri\ica\res_frontal25\'; 
load(strcat(mypath,'labels.mat'));
ic_list=[15,25,18,14,21,16,5,22];

% list of task sessions
ses_list={'simon20','simon50','stroop20','stroop50'};

%%% set parameters
% define epoch for hemodynamic response
t = 0:0.1:9.9;
% define matrix for hdr & its peak
hdr = nan(height(T),length(t),length(ic_list));
peak = nan(height(T),length(ic_list));
% time range for peak extraction
peakTime = t > 3 & t < 5;

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
            hdrSes(tr,:,:)=itc(tr_list(tr):tr_list(tr)+length(t)-1,ic_list) - itc(tr_list(tr),ic_list);
        end
        
        hdr(SesIdx,:,:) = hdrSes;
        peak(SesIdx,:) = mean(hdrSes(:,peakTime,:),2);
         
    end

end


%% 
Tmri = [T array2table(peak,'VariableNames',labels(ic_list))];

col1 = find(strcmp(Tmri.Properties.VariableNames,'SMA'));
col_stimRT = find(strcmp(Tmri.Properties.VariableNames,'stimRT'));
col_normRT = find(strcmp(Tmri.Properties.VariableNames,'normRT'));

% change here between 'conflict' and 'accuracy' conditions
trials = ( T.type=='incon') & T.stimACC==1  & T.badep==1 ; 
% trials = T.type=='incon'  & T.stimRT~=0 & T.badep==1; T.type=='con' |

Tmri = Tmri(trials,:);
Tmri.type   = removecats(Tmri.type);
Tmri.freq   = removecats(Tmri.freq); 
Tmri.task   = removecats(Tmri.task); 
Tmri.isicat = removecats(Tmri.isicat);

% zscore for each subject
for sb=sb_in    
    
   for col = col1:width(Tmri)
       Tmri{Tmri.Subject==sb,col} = zscore(Tmri{Tmri.Subject==sb,col});
   end

   Tmri{Tmri.Subject==sb,col_normRT} = zscore(Tmri.stimRT(Tmri.Subject==sb));
    
end
%%
% calculate correlation maps - time consuming
tic
%'SMA','area8',,'area10','AIC','PIC'
ic_names = {'preSMA','ACC','area9'};
ltyp = 2;

estim = nan(num_frex,length(time(ltyp).t),length(ic_names)); %
pval = nan(num_frex,length(time(ltyp).t),length(ic_names)); %


for ic = 1:length(ic_names)
for f = 1:num_frex
    parfor tp = 1:length(time(ltyp).t)
        psd = nan(height(T),1); 
        for sb = sb_in
            psd(T.Subject==sb,1) = tf(sb,ltyp).data(f,tp,:);
            psd(T.Subject==sb,1) = zscore(psd(T.Subject==sb,1));
        end
        psd = psd(trials ,:);
        Tall = [Tmri table(psd)];

        lme = fitlme(Tall,strcat('psd ~ ', ic_names{ic} ,'+ (1|Block) + (1|isicat) + (1|Subject)'));
        
        estim(f,tp,ic)= lme.Coefficients(2,2).Estimate;
        pval(f,tp,ic) = lme.Coefficients(2,6).pValue;
    end
end
end
toc
%%
hFig = figure(1); clf
set( hFig, 'units','normalized','outerposition',[0.3 0.5 0.3 0.5]);

timlim = [time(ltyp).t(1) time(ltyp).t(end)];

for es = 1:length(ic_names)
    subplot(2,2,es)
    
    pcor = fdr_bh(squeeze(pval(:,:,es)),0.1); 
    
    contourf(time(ltyp).t,frex,squeeze(estim(:,:,es)),50,'linecolor','none');  
    hold on
    contour(time(ltyp).t,frex,pcor,1,'LineWidth',3,'linecolor','k');

    xlabel('time [ms]'); ylabel('frequency [Hz]'); title(ic_names(es))
    set(gca,'clim',[-.05 .05],'xlim',timlim,'ydir','norm'); colorbar; colormap(jet(256))
end
