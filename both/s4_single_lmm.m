% Script written 1.02.2021 by Ewa Beldzik 
% traditional EEG and fMRI analyses
% coupling EEG & fMRI using 'single' LME

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
time(1).t = -300:10:1200;
timeStimIdx = dsearchn(EEG.times',time(1).t'); 

% response-locked time epoch
time(2).t = -800:10:400;

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
    tf(sb,2).data = nan(num_frex,length(time(2).t),ntrials);   
    
    for tr=find(typIdx)'
        timeBresp = time(2).t + Tsb.stimRT(tr);
        timeBrespIdx = dsearchn(EEG.times',timeBresp');
        tf(sb,2).data(:,:,tr) = tf_temp(:,timeBrespIdx,tr);
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

%% EEG export mean frequency bands from tf maps
band_list = {'delta','theta','alpha','beta'};
f = {frex <= 3, frex >= 4 & frex <= 8, frex >= 9 & frex <= 12, frex >= 13 & frex <= 30};
timeBresp = -300:0; timeBrespIdx = unique(dsearchn(time(2).t',timeBresp'));
timeAresp = 0:300; timeArespIdx = unique(dsearchn(time(2).t',timeAresp'));

% create new table with EEG data
Teeg = [];

% initiate results
psdS = []; 
thetaMax    = nan(height(T),1); 
thetaBefor  = nan(height(T),1); 
thetaAfter  = nan(height(T),1); 

for band = 1:length(f)
    psdS(band).data = nan(height(T),1);
    
    for sb = sb_in 
        
          % extract stimulus-locked mean
          psdS(band).data(T.Subject==sb,1) = mean(mean(tf(sb,1).data(f{band},time(1).t > 0 ,:)));          
          
          if strcmp(band_list{band},'theta')
              % extract max theta, resp-locked theta before & after resp
                thetaMax(T.Subject==sb,1) = max(mean(tf(sb,1).data(f{band},time(1).t > 0 ,:)));               
                thetaBefor(T.Subject==sb,1) = mean(mean(tf(sb,2).data(f{band},timeBrespIdx,:))); 
                thetaAfter(T.Subject==sb,1) = mean(mean(tf(sb,2).data(f{band},timeArespIdx,:)));  
          end           
          
    end 
    Teeg = [Teeg array2table(psdS(band).data,'VariableNames',strcat(band_list(band),'S'))]; 

end
Teeg = [Teeg table(thetaMax) table(thetaBefor) table(thetaAfter)];

%% make table   
Tall = [T array2table(peak,'VariableNames',labels(ic_list)) Teeg]; 

% change here between 'conflict' and 'accuracy' conditions
% Tall = Tall((Tall.type=='con' | Tall.type=='incon') & Tall.accur=='hit' & Tall.badep==1,:); 
Tall = Tall(Tall.type=='incon' & Tall.stimRT~=0 & Tall.badep==1 ,:);

Tall.type   = removecats(Tall.type);
Tall.freq   = removecats(Tall.freq); 
Tall.task   = removecats(Tall.task); 
Tall.isicat = removecats(Tall.isicat);
Tall.accur = removecats(Tall.accur);
col1        = find(strcmp(Tall.Properties.VariableNames,'SMA'));
col_stimRT  = find(strcmp(Tall.Properties.VariableNames,'stimRT'));
col_normRT  = find(strcmp(Tall.Properties.VariableNames,'normRT'));
col_bl      = find(strcmp(Tall.Properties.VariableNames,'Block'));
col_tr      = find(strcmp(Tall.Properties.VariableNames,'Trial'));

% zscore for each subject
for sb=sb_in    
    
    for col = col1:width(Tall)
         Tall{Tall.Subject==sb,col} = zscore(Tall{Tall.Subject==sb,col});
    end

    Tall{Tall.Subject==sb,col_normRT} = zscore(Tall.stimRT(Tall.Subject==sb));
    Tall{Tall.Subject==sb,col_bl}     = zscore(Tall.Block(Tall.Subject==sb));
    Tall{Tall.Subject==sb,col_tr}     = zscore(Tall.Trial(Tall.Subject==sb));
end

%% Behaioural model
hFig = figure(1); clf; hold on
set( hFig, 'units','normalized','outerposition',[1 0.66 0.22 0.3]);

lme = fitlme(Tall,'stimRT ~ type*freq + task + dur + (1|Subject)'); 

% set order for plotting 
sortEs = [3,2,6,4,5];
names = lme.CoefficientNames(sortEs);
estim = lme.Coefficients(sortEs,2).Estimate.*[1;-1;1;1;1];
upCI  = lme.Coefficients(sortEs,8).Upper.*[1;-1;1;1;1];
lowCI = lme.Coefficients(sortEs,7).Lower.*[1;-1;1;1;1];
wid = 0.3;

for i = 1:length(estim)     
    if strcmp(names(i),'task_simon')
        rectangle('Position',[i-wid,upCI(i),2*wid,lowCI(i)-upCI(i)], 'FaceColor',[.85 .7 .85])
    else
        rectangle('Position',[i-wid,lowCI(i),2*wid,upCI(i)-lowCI(i)], 'FaceColor',[.85 .7 .85])

    end
    
    plot([i-wid i+wid],[estim(i) estim(i)],'color',[.5  .1  .4],'linewidth',2);     
    
end
plot([0 length(estim)+0.5],[0 0],'color','k','linewidth',1);
xticks(1:length(estim)); xlim([0.5 length(estim)+0.5]); xticklabels(names)
ylim([-25 80]);

%% average and plot fMRI HDR
hFig = figure(2); clf; hold on;
set( hFig, 'units','normalized','outerposition',[0.7 0 0.1 1]);

hdrAve=[]; hdrTyp = []; errBar=[];
types = {'con','incon'}; colors = {'c','b'};
% types = {'hit','err'}; colors = {'b','r'};

% average HDR for each subject and conflict type
for trialtyp=1:length(types) 
    for sb=sb_in    
        % change here between 'conflict' and 'accuracy' conditions
        typIdx = T.Subject==sb & T.type==types(trialtyp) & T.stimACC==1 & T.badep==1;
%       typIdx = T.Subject==sb & T.accur==types(trialtyp) & T.type=='incon' & T.stimRT~=0 & T.badep==1;
        hdrTyp(trialtyp).data(sb,:,:) = mean(hdr(typIdx,:,:));
    end
end
pval = nan(length(t),1);
corIdx = t > 1;

for ic=1:length(ic_list)
    
    subplot(8,1,ic)
    [~,pval(corIdx)] = ttest(hdrTyp(1).data(:,corIdx,ic),hdrTyp(2).data(:,corIdx,ic));  
    pcor = 100*fdr_bh(pval,0.01); pcor(pcor==0)=NaN;  
    
    area(t,pcor,-100,'LineStyle','none','FaceColor',[.87 .87 .87]);
    hold on
    for trialtyp=1:length(types)
        hdrAve(trialtyp,:) = nanmean(hdrTyp(trialtyp).data(sb_in,:,ic));
        errBar(trialtyp,:) = nanstd(hdrTyp(trialtyp).data(sb_in,:,ic))/sqrt(length(sb_in));     
        shadedErrorBar(t,hdrAve(trialtyp,:),errBar(trialtyp,:),'lineprops',colors{trialtyp}); 
    end
    
    ylim([-0.25 0.25]); %*4
    title(strcat('ic',num2str(ic_list(ic)),'-',labels(ic_list(ic))))

end
xlabel('time [s]');

%% average and plot TF maps
figure(3); clf
% change here between 'conflict' and 'accuracy' conditions
types = {'con','incon'};
% types = {'hit','err'};

for locktyp=1:2

    % averaged TR maps
    tfTyp = []; tfAll = [];
    for sb=sb_in  
        Tsb=T(T.Subject==sb,:); 
    
        %for all trials
        typIdx=(Tsb.type=='con' | Tsb.type=='incon') & Tsb.stimACC==1 & Tsb.badep == 1;
        tfAll(sb,:,:) = mean(tf(sb,locktyp).data(:,:,typIdx),3);

        %for each trials type           
        for trialtyp=1:length(types)  
            % change here between 'conflict' and 'accuracy' conditions
            typIdx = Tsb.type==types(trialtyp) & Tsb.stimACC==1 & Tsb.badep == 1;
            % typIdx = Tsb.accur==types(trialtyp) & Tsb.type=='incon' & Tsb.stimRT~=0 & Tsb.badep == 1;
            tfTyp(trialtyp).data(sb,:,:) = mean(tf(sb,locktyp).data(:,:,typIdx),3); 
        end  
      
    end
    
    % perform t tests
    p = nan(num_frex,length(time(locktyp).t)); pcor =  nan(size(p)); 
    for fi = 1:num_frex
        [~, p(fi,:)] = ttest(tfTyp(1).data(sb_in,fi,:),tfTyp(2).data(sb_in,fi,:));
    end
    pcor = fdr_bh(p,0.001);
    
    % plot results    
    % topoplot(mean(ic_topo),EEG.chanlocs);  
    timlim = [time(locktyp).t(1) time(locktyp).t(end)];
    
    subplot(2,2,2*locktyp-1)
    contourf(time(locktyp).t,frex,squeeze(mean(tfAll(sb_in,:,:))),50,'linecolor','none');   
    
    xlabel('time [ms]'); ylabel('frequency [Hz]'); title('all correct')
    set(gca,'clim',[-.08 .08],'xlim',timlim,'ydir','norm'); xticks(-600:300:1200); colormap(jet(256)); % colorbar; 

    subplot(2,2,2*locktyp)
    diffmap = squeeze(mean(tfTyp(2).data(sb_in,:,:))) - squeeze(mean(tfTyp(1).data(sb_in,:,:)));
    maxVal = max(max(diffmap));    
    contourf(time(locktyp).t,frex,diffmap,50,'linecolor','none'); 
    hold on
     contour(time(locktyp).t,frex,pcor,1);
%     window = zeros(size(p)); window(frex >= 4 & frex <= 8,time(2).t > 0 & time(2).t < 300)=1; 
%     contour(time(locktyp).t,frex,window,1,'LineStyle','--');
    
    xlabel('time [ms]'); ylabel('frequency [Hz]'); title(strcat(types(2),{' - '},types(1)))
    set(gca,'clim',[-.03 .03],'xlim',timlim,'ydir','normal'); xticks(-600:300:1200); colormap(jet(256)); % *6
end

%% EEG-fMRI single model
hFig = figure(5); clf
set( hFig, 'units','normalized','outerposition',[0.1 0.3 0.5 0.7]);

ExtType = {'thetaS' 'thetaMax' 'thetaBefor' 'thetaAfter'};

wid = 0.3;
for tt=1:length(ExtType)
    
    %%% run LMM with behavioral variables
    % change here between 'conflict' and 'accuracy' conditions
    %  lme = fitlme(Tall, strcat(ExtType{tt}, '~ normRT + type + (1|Block) + (1|isicat) + (1|Subject)')); % 
    lme = fitlme(Tall, strcat(ExtType{tt}, '~ normRT + accur + (1|Block) + (1|isicat) + (1|Subject)'));
    
    % read LMM results     
    names = lme.CoefficientNames(2:3);
    pval  = lme.Coefficients(2:3,6).pValue;
    estim = lme.Coefficients(2:3,2).Estimate;
    upCI  = lme.Coefficients(2:3,8).Upper;
    lowCI = lme.Coefficients(2:3,7).Lower;    
    
    subplot(4,7,7*tt-4);   
    plot([0.3 length(estim)+0.7],[0 0],'color','k','linewidth',1);
    for i = 1:length(estim)
        % indicate color of the bars
        if pval(i) < 0.05
            colorCI = [.85 .7 .85]; %purple
            colorEs = [.5  .1  .4];
        else    
            colorCI = [.8 .8 .8]; %grey
            colorEs = [.4 .4 .4];
        end

        hold on
        rectangle('Position',[i-wid,lowCI(i),2*wid,upCI(i)-lowCI(i)],'FaceColor',colorCI)
        plot([i-wid i+wid],[estim(i) estim(i)],'color',colorEs,'linewidth',3); 

    end    
    xlim([0.5 length(estim)+0.5]); xticks(1:length(estim)); xticklabels(names);
    ylim([0 .1])
    
    
    %%% run LMM for brain networks
    lme = fitlme(Tall, strcat(ExtType{tt},... 
    '~ SMA + preSMA + ACC + area8 + area9 + area10 + AIC + PIC + deltaS + alphaS + betaS + (1|Block) + (1|isicat) + (1|Subject)'));

    % read LMM results    
    EsIdx = 2:lme.NumCoefficients -3;
    names = lme.CoefficientNames(EsIdx);
    pval  = lme.Coefficients(EsIdx,6).pValue;2*i-1;
    estim = lme.Coefficients(EsIdx,2).Estimate;
    upCI  = lme.Coefficients(EsIdx,8).Upper;
    lowCI = lme.Coefficients(EsIdx,7).Lower;
    
    subplot(4,2,2*tt);
    plot([0.3 length(estim)+0.7],[0 0],'color','k','linewidth',1);
    for i = 1:length(estim)
        
        % indicate color of the bars
        if pval(i) < 0.05 && estim(i) > 0
            colorCI = [.9 .3 .3]; %red
            colorEs = [.5  0  0];
        elseif pval(i) < 0.05 && estim(i) < 0
            colorCI = [.3 .6 .8]; %blue
            colorEs = [ 0 .3 .5];
        else    
            colorCI = [.8 .8 .8]; %grey
            colorEs = [.4 .4 .4];
        end
        
        hold on
        rectangle('Position',[i-wid,lowCI(i),2*wid,upCI(i)-lowCI(i)],'FaceColor',colorCI)
        plot([i-wid i+wid],[estim(i) estim(i)],'color',colorEs,'linewidth',3); 

    end
    xlim([0.5 length(estim)+0.5]); xticks(1:length(estim)); xticklabels(names);
    ylim([-.05 .05]); yticks(-0.05:0.05:0.05);
end

%% Additional analyses - block and trial effects
hFig = figure(6); clf
set( hFig, 'units','normalized','outerposition',[0.1 0.3 0.3 0.4]);

var2plot = {'thetaS','area9'};
wid = 0.28;
pval =[];

% run LMM for stim RT
fitlme(Tall, 'stimRT ~ type + Block + Trial + (1|isicat) + (1|Subject)')

% run LMM for brain variables
for var = 1:length(var2plot)
    lme = fitlme(Tall, strcat(var2plot{var},' ~ normRT + type + Block + Trial + (1|isicat) + (1|Subject)'));
    
    % read LMM results     
    names = lme.CoefficientNames([5,4,3,2]);
    pval(:,var)  = lme.Coefficients([5,4,3,2],6).pValue;
    estim = lme.Coefficients([5,4,3,2],2).Estimate;
    upCI  = lme.Coefficients([5,4,3,2],8).Upper;
    lowCI = lme.Coefficients([5,4,3,2],7).Lower;    
    
    subplot(1,2,var);   
    plot([0.3 length(estim)+0.7],[0 0],'color','k','linewidth',1);

    for i = 1:length(estim)
        % indicate color of the bars
        if i < 3
            colorCI = [.85 .7 .85]; % purple
            colorEs = [.5  .1  .4];
        else    
            colorCI = [1 .8 0]; % orange
            colorEs = [.8 .5 0];
        end
    
        hold on
        rectangle('Position',[i-wid,lowCI(i),2*wid,upCI(i)-lowCI(i)],'FaceColor',colorCI)
        plot([i-wid i+wid],[estim(i) estim(i)],'color',colorEs,'linewidth',3); 
    
    end    
    xlim([0.5 length(estim)+0.5]); xticks(1:length(estim)); xticklabels(names);
    ylim([-.08 .08])
    set(gca,'fontsize', 13)
end
