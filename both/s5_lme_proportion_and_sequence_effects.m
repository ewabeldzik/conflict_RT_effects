% Script written 25.07.2022 by Ewa Beldzik 
% conducts the LME analyses aimed to verify the RT, proportion congruency,
% and congruency sequence effects in EEG and fMRI 'conflict markers'
% 
% Results are presented in the Figure 4 in the manuscript entitled:
% "A thin line between conflict and reaction time effects on EEG and fMRI
% brain signals" (2023) by Ewa Beldzik & Markus Ullsperger


clear; clc

load('E:\SEF2\beh\beh_sef.mat');

a = T.before=='blank';
b = a;
b(1:2)=[];
b(end+2) = false;
T.before(a) = T.type(b);
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

% response-locked time epoch
timeEEG = -800:10:400;

%%% wavelet parameters
num_frex = 15;
range_frex = [ 4 8 ];
range_cycl = [ 3 4 ];
%logarithmically spaced frequency
frex  = logspace(log10(range_frex(1)),log10(range_frex(2)),num_frex);
nCycs = logspace(log10(range_cycl(1)),log10(range_cycl(2)),num_frex); 

% initialize output time-frequency data
tf = [];
thetapow = nan(height(T),length(timeEEG));


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
    
    % prepare resp-locked tf maps
    Tsb = T(T.Subject==sb,:);  
    typIdx = Tsb.stimRT~=0;     
    tf(sb).data = nan(num_frex,length(timeEEG),ntrials);   
    
    for tr=find(typIdx)'
        timeBresp = timeEEG + Tsb.stimRT(tr);
        timeBrespIdx = dsearchn(EEG.times',timeBresp');
        tf(sb).data(:,:,tr) = tf_temp(:,timeBrespIdx,tr);
    end
    
    thetapow(T.Subject==sb,:) = squeeze(mean(tf(sb).data))';
end

%% fMRI analysis
% load fMRI data
mypath  = 'E:\SEF2\mri\ica\res_frontal25\'; 
load(strcat(mypath,'labels.mat'));
ic_list=[18,25];

% list of task sessions
ses_list={'simon20','simon50','stroop20','stroop50'};

%%% set parameters
% define epoch for hemodynamic response
timeBOLD = 0:0.1:9.9;
% define matrix for hdr & its peak
hdr = nan(height(T),length(timeBOLD),length(ic_list));


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
        hdrSes=nan(length(tr_list),length(timeBOLD),length(ic_list));
        
        %extract HDR of each trial
        for tr=1:length(tr_list)  
            hdrSes(tr,:,:) = itc(tr_list(tr):tr_list(tr)+length(timeBOLD)-1,ic_list) - itc(tr_list(tr),ic_list);
        end
        
        hdr(SesIdx,:,:) = hdrSes;
         
    end

end

%% Prepare variable 'brain' with all brain measures of interest

brainVars = {'theta' 'ACC' 'preSMA'};
brain = [];

for sb = sb_in     
    
    brain(1).tc(sb,:,:) = squeeze(thetapow(T.Subject==sb,:))';    
    brain(2).tc(sb,:,:) = squeeze(hdr(T.Subject==sb,:,find(ismember(labels(ic_list),'ACC'))))';    
    brain(3).tc(sb,:,:) = squeeze(hdr(T.Subject==sb,:,find(ismember(labels(ic_list),'preSMA'))))';

end


%% Run LME analyses for RT effect
for bvar = 1:length(brainVars) 

    if bvar == 1; time2plot = timeEEG; else time2plot = timeBOLD; end
    for tp = 2:length(time2plot)

        % extract brain value at each time point (tp)
        brainVal = nan(height(T),1);
        for sb = sb_in
            brainVal(T.Subject==sb)  = squeeze(brain(bvar).tc(sb,tp,:));    
        end
    
        % create new table with correct trials only
        Tall = [T array2table(brainVal,'VariableNames',brainVars(bvar))]; 
        Tx = Tall((Tall.type=='con' | Tall.type=='incon') & (Tall.before=='con' | Tall.before=='incon') & ...
                   Tall.accur=='hit'  & Tall.badep==1 & Tall.Trial>1 & Tall.stimRT~=0,:); 
        
        Tx.type   = removecats(Tx.type);
        Tx.isicat = removecats(Tx.isicat);
        col1      = find(strcmp(Tx.Properties.VariableNames,'theta'));
        
        % zscore for each subject
        for sb=sb_in            
            Tx.normRT(Tx.Subject==sb) = zscore(Tx.stimRT(Tx.Subject==sb));
            for col = col1:col1 + width(Tall) - width(T) -1
                 Tx{Tx.Subject==sb,col} = zscore(Tx{Tx.Subject==sb,col});
            end   
        end
    
        %%% run LMM with behavioral variables
        lme = fitlme(Tx, strcat(brainVars{bvar}, '~ type + normRT + (1 + type + normRT |Subject)')); %
    
        % read LMM results     
        brain(bvar).pval(:,tp)  = lme.Coefficients(:,6).pValue;
        brain(bvar).estim(:,tp) = lme.Coefficients(:,2).Estimate;
        brain(bvar).upCI(:,tp)  = lme.Coefficients(:,8).Upper;
        brain(bvar).lowCI(:,tp) = lme.Coefficients(:,7).Lower;
        

    end
    
end
[~,ConIdx] = ismember('type_incon',lme.CoefficientNames);


%% plot LME analyses for RT effect
hFig = figure(1); clf
set( hFig, 'units','normalized','outerposition',[.1 .1 0.4 0.8]);

% define colors
VarList = {'type_incon','freq_congr20','type_incon:freq_congr20','before_incon','type_incon:before_incon','normRT'};

ColorListCI = {[.49 .73 .88];...  % blue
               [.72 .7  .86];...  % purple 
               [.94 .68 .82];...  % pink
               [.55 .82 .88];...  % turquois
               [.65 .91 .9];...   % green
               [.9  .71 .3]};     % yellow


ColorListEs = {[.15 .42 .62];...  % blue
               [.5  .4  .6];...   % purple 
               [.84 .42 .73];...  % pink
               [.1  .55 .61];...  % turquois
               [.12 .66 .60];...  % green
               [.6  .35 .2]};     % yellow


for bvar = 1:length(brainVars) 
    errBar = brain(bvar).upCI - brain(bvar).estim;

    subplot(4,3,3*bvar+1)
    if bvar == 1; time2plot = timeEEG; else time2plot = timeBOLD; end

    for i = 2:3    
        hold on
        shadedErrorBar(time2plot,brain(bvar).estim(i,:),errBar(i,:),'lineProps',{'markerfacecolor',[1 1 1]});
        xlim('tight'); 
    end
    % changing the color of RT var (weird but works)
    shadedErrorBar(time2plot,brain(bvar).estim(i,:),errBar(i,:),'lineProps',{'markerfacecolor',[1 1 1]});
    shadedErrorBar(time2plot,brain(bvar).estim(i,:),errBar(i,:),'lineProps',{'markerfacecolor',[1 1 1]});
    plot(time2plot,brain(bvar).estim(i,:),'color',[.6  .35 .2]);

    if bvar==1;  xticks(-600:300:300); end
end
%legend(names(2:3))

%% create table with peak 'congruency' condition

% find max 'congruency effect' and extract brain values
brainVal = nan(height(T),length(brainVars));
for bvar = 1:length(brainVars) 
    
    [~, brain(bvar).max] = max(brain(bvar).estim(ConIdx,:));
    for sb = sb_in
        brainVal(T.Subject==sb,bvar)  = squeeze(mean(brain(bvar).tc(sb,brain(bvar).max-5:brain(bvar).max+5,:)));    
    end
end

% create new table  with 'previous' trial (for Gratton analysis)
Tall = [T array2table(brainVal,'VariableNames',brainVars)]; 

% select correct trials
Tx = Tall((Tall.type=='con' | Tall.type=='incon') & (Tall.before=='con' | Tall.before=='incon') & ...
           Tall.accur=='hit'  & Tall.badep==1 & Tall.Trial>1 & Tall.stimRT~=0,:);

Tx.type   = removecats(Tx.type);
Tx.before  = removecats(Tx.before);
Tx.freq   = removecats(Tx.freq);
Tx.isicat = removecats(Tx.isicat);
col1      = find(strcmp(Tx.Properties.VariableNames,'theta'));

% zscore for each subject
for sb=sb_in            
    Tx.normRT(Tx.Subject==sb) = zscore(Tx.stimRT(Tx.Subject==sb));
    for col = col1:col1 + length(brainVars) -1
         Tx{Tx.Subject==sb,col} = zscore(Tx{Tx.Subject==sb,col});
    end   
end


%% plot the 2nd column - proportion congruency effect
wid = 0.25;

%%% run LMM on RTs
lme = fitlme(Tx, 'stimRT ~  freq*type +  (1 + freq*type |Subject)'); %
FixEff = size(lme.Coefficients,1);

% read LMM results     
names = lme.CoefficientNames(2:FixEff);
pval  = lme.Coefficients(2:FixEff,6).pValue;
estim = lme.Coefficients(2:FixEff,2).Estimate;
upCI  = lme.Coefficients(2:FixEff,8).Upper;
lowCI = lme.Coefficients(2:FixEff,7).Lower;

subplot(4,3,2)
plot([0.3 length(estim)+1.7],[0 0],'color','k','linewidth',1);

for i = 1:length(estim)

        % indicate color of the bars
        if pval(i) < 0.05
            colorCI = ColorListCI{find(strcmp(VarList,names(i)))}; 
            colorEs = ColorListEs{find(strcmp(VarList,names(i)))}; 
        else       
            colorCI = [.8 .8 .8]; %grey
            colorEs = [.4 .4 .4];
        end

    hold on
    rectangle('Position',[i-wid,lowCI(i),2*wid,upCI(i)-lowCI(i)],'FaceColor',colorCI)
    plot([i-wid i+wid],[estim(i) estim(i)],'color',colorEs,'linewidth',3); 

end    

xlim([0.5 length(estim)+1.5]); xticks(1:length(estim)); %xticklabels(names);
ylim([min(lowCI)*1.5 max(upCI)*1.1])


%%% run LMM on brain measures
for bvar = 1:length(brainVars)

    
    lme = fitlme(Tx, strcat(brainVars{bvar}, '~  freq*type + normRT + (1 + freq*type + normRT|Subject)')); %
    FixEff = size(lme.Coefficients,1);

    % read LMM results     
    names = lme.CoefficientNames([2,3,5,4]);
    pval  = lme.Coefficients([2,3,5,4],6).pValue;
    estim = lme.Coefficients([2,3,5,4],2).Estimate;
    upCI  = lme.Coefficients([2,3,5,4],8).Upper;
    lowCI = lme.Coefficients([2,3,5,4],7).Lower;

    subplot(4,3,3*bvar+2)
    plot([0.3 length(estim)+0.7],[0 0],'color','k','linewidth',1);

    for i = 1:length(estim)

        % indicate color of the bars
        if pval(i) < 0.05
            colorCI = ColorListCI{find(strcmp(VarList,names(i)))}; 
            colorEs = ColorListEs{find(strcmp(VarList,names(i)))}; 
        else       
            colorCI = [.8 .8 .8]; %grey
            colorEs = [.4 .4 .4];
        end

        hold on
        rectangle('Position',[i-wid,lowCI(i),2*wid,upCI(i)-lowCI(i)],'FaceColor',colorCI)
        plot([i-wid i+wid],[estim(i) estim(i)],'color',colorEs,'linewidth',3); 

    end    

    xlim([0.5 length(estim)+0.5]); xticks(1:length(estim)); %xticklabels(names);
    ylim([min(lowCI)*1.3 max(upCI)*1.1])

end

%% plot the 3rd column - congruency sequence effect

%%% run LMM on RTs
lme = fitlme(Tx, 'stimRT ~ before*type +  (1 + before*type |Subject)'); %
FixEff = size(lme.Coefficients,1);

% read LMM results     
names = lme.CoefficientNames(2:FixEff);
pval  = lme.Coefficients(2:FixEff,6).pValue;
estim = lme.Coefficients(2:FixEff,2).Estimate;
upCI  = lme.Coefficients(2:FixEff,8).Upper;
lowCI = lme.Coefficients(2:FixEff,7).Lower;

subplot(4,3,3)
plot([0.3 length(estim)+0.7],[0 0],'color','k','linewidth',1);

for i = 1:length(estim)

    % indicate color of the bars
    if pval(i) < 0.05
        colorCI = ColorListCI{find(strcmp(VarList,names(i)))}; 
        colorEs = ColorListEs{find(strcmp(VarList,names(i)))}; 
    else    
        colorCI = [.8 .8 .8]; %grey
        colorEs = [.4 .4 .4];
    end

    hold on
    rectangle('Position',[i-wid,lowCI(i),2*wid,upCI(i)-lowCI(i)],'FaceColor',colorCI)
    plot([i-wid i+wid],[estim(i) estim(i)],'color',colorEs,'linewidth',2); 

end    
xlim([0.5 length(estim)+1.5]); xticks(1:length(estim)); %xticklabels(names);
ylim([min(lowCI)*2 max(upCI)*1.1])


%%% run LMM with behavioral variables
for bvar=1:length(brainVars)

    
    lme = fitlme(Tx, strcat(brainVars{bvar}, '~ before*type + normRT +  (1 + before*type + normRT |Subject)')); %
    FixEff = size(lme.Coefficients,1);

    % read LMM results     
    names = lme.CoefficientNames([2,3,5,4]);
    pval  = lme.Coefficients([2,3,5,4],6).pValue;
    estim = lme.Coefficients([2,3,5,4],2).Estimate;
    upCI  = lme.Coefficients([2,3,5,4],8).Upper;
    lowCI = lme.Coefficients([2,3,5,4],7).Lower;

    subplot(4,3,3*bvar+3)
    plot([0.3 length(estim)+0.7],[0 0],'color','k','linewidth',1);

    for i = 1:length(estim)

        % indicate color of the bars
        if pval(i) < 0.05
            colorCI = ColorListCI{find(strcmp(VarList,names(i)))}; 
            colorEs = ColorListEs{find(strcmp(VarList,names(i)))}; 
        else       
            colorCI = [.8 .8 .8]; %grey
            colorEs = [.4 .4 .4];
        end

        hold on
        rectangle('Position',[i-wid,lowCI(i),2*wid,upCI(i)-lowCI(i)],'FaceColor',colorCI)
        plot([i-wid i+wid],[estim(i) estim(i)],'color',colorEs,'linewidth',2); 

    end    

    xlim([0.5 length(estim)+0.5]); xticks(1:length(estim)); %xticklabels(names);
    ylim([min(lowCI)*1.3 max(upCI)*1.1])

end

