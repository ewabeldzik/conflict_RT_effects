% Script written 25.07.2022 by Ewa Beldzik 
% conducts the LME analyses aimed to verify the RT, proportion congruency,
% and congruency sequence effects in EEG and fMRI 'conflict markers'
% 
% Results are presented in the Figure 4 in the manuscript entitled:
% "A thin line between conflict and reaction time effects on EEG and fMRI
% brain signals" (2023) by Ewa Beldzik & Markus Ullsperger


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
con20 = T.freq == 'congr20' & T.type=='con';
con50 = T.freq == 'congr50' & T.type=='con';
incon50 = T.freq == 'congr50' & T.type=='incon';
incon20 = T.freq == 'congr20' & T.type=='incon';
cc = T.before == 'con' & T.type=='con';
ic = T.before == 'incon' & T.type=='con';
ci = T.before == 'con' & T.type=='incon';
ii = T.before == 'incon' & T.type=='incon';
T = [T table(con50) table(con20)  table(incon50) table(incon20) table(ic) table(cc) table(ii) table(ci) ];

% load fMRI data
mypath  = 'E:\SEF2\mri\ica\res_frontal25\'; 
load(strcat(mypath,'labels.mat'));
ic_list = [18, 25];

% list of task sessions
ses_list = {'simon20','simon50','stroop20','stroop50'};

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
            hdrSes(tr,:,:) = itc(tr_list(tr):tr_list(tr)+length(timeBOLD)-1,ic_list);
        end
        
        hdr(SesIdx,:,:) = hdrSes;
         
    end

end



%% z-score values
trSelIdx = (T.type=='con' | T.type=='incon') & (T.before=='con' | T.before=='incon') & ...
                   T.accur=='hit'  & T.badep==1 & T.Trial>1 & T.stimRT~=0; 

Tx = T(trSelIdx,:); 
Tx.type   = removecats(Tx.type);
Tx.before  = removecats(Tx.before);
Tx.freq   = removecats(Tx.freq);

% zscore for each subject
for sb=sb_in            
    Tx.normRT(Tx.Subject==sb) = zscore(Tx.stimRT(Tx.Subject==sb));
end



%% Prepare variable 'brain' with all brain measures of interest

hdr2=hdr(trSelIdx,:,:);
for bvar = 1:2
    for sb=sb_in
        sbIdx = Tx.Subject==sb;
        hdr2(sbIdx,:,bvar) = zscore(hdr2(sbIdx,:,bvar));
    end
hdr2(:,:,bvar) = hdr2(:,:,bvar) - hdr2(:,1,bvar);
end

%% Run LME analyses
figure(1); clf
colorFreqIdx  = {[189 21 29]/256, [126 47 142]/256, [71 28 11]/256,  [168 148 22]/256};
colorBeforIdx = {[77 87 199]/256, [0 114 189]/256,  [34 209 39]/256, [10 143 79]/256 };


brainVars = {'ACC' 'preSMA'};
brain = [];
time2plot = timeBOLD;

bvar = 1;
for tp = 2:length(time2plot)

    % extract brain value at each time point (tp)
    brainVal = nan(height(Tx),1);
    for sb = sb_in
        brainVal(Tx.Subject==sb)  = squeeze(hdr2(Tx.Subject==sb,tp,bvar));    
    end

    % create new table with correct trials only
    Tz = [Tx array2table(brainVal,'VariableNames',brainVars(bvar))]; 

    %%% run LMM with behavioral variables
    % lme = fitlme(Tz, strcat(brainVars{bvar}, '~ -1 + con50 + con20 + incon50 + incon20 + normRT + (-1 + con20 + con50 + incon20 + incon50 + normRT |Subject)')); %
    lme = fitlme(Tz, strcat(brainVars{bvar}, '~ -1 + con20 + con50 + incon20 + incon50 + (-1 + con20 + con50 + incon20 + incon50 |Subject)')); %')); %

    % read LMM results     
    brain(bvar).pval(:,tp)  = lme.Coefficients(:,6).pValue;
    brain(bvar).estim(:,tp) = lme.Coefficients(:,2).Estimate;
    brain(bvar).upCI(:,tp)  = lme.Coefficients(:,8).Upper;
    brain(bvar).lowCI(:,tp) = lme.Coefficients(:,7).Lower;
    

end

%%% plot ACC
subplot(1,2,bvar)
cIdx=1;
[~,ConIdx] = ismember('con50_1',lme.CoefficientNames);
errBar = brain(bvar).upCI - brain(bvar).estim;
for i = ConIdx:ConIdx+3   
    hold on
    shadedErrorBar(time2plot,brain(bvar).estim(i,:),errBar(i,:),'lineProps',{'-','color',colorFreqIdx{cIdx}, 'LineWidth',1.5});
    xlim('tight'); 
    cIdx = cIdx +1;
end

legend(lme.CoefficientNames(ConIdx:ConIdx+3))

%%
bvar = 2;
for tp = 2:length(time2plot)

    % extract brain value at each time point (tp)
    brainVal = nan(height(Tx),1);
    for sb = sb_in
        brainVal(Tx.Subject==sb)  = squeeze(hdr2(Tx.Subject==sb,tp,bvar));    
    end

    % create new table with correct trials only
    Tz = [Tx array2table(brainVal,'VariableNames',brainVars(bvar))]; 

    %%% run LMM with behavioral variables
    lme2 = fitlme(Tz, strcat(brainVars{bvar}, '~ -1 + cc + ci + ic + ii + normRT + (-1 + cc + ci + ic + ii + normRT |Subject)')); %
     % lme2 = fitlme(Tz, strcat(brainVars{bvar}, '~ -1 + cc + ci + ic + ii + (-1 + cc + ci + ic + ii |Subject)')); %')); % 

    % read LMM results     
    brain(bvar).pval(:,tp)  = lme2.Coefficients(:,6).pValue;
    brain(bvar).estim(:,tp) = lme2.Coefficients(:,2).Estimate;
    brain(bvar).upCI(:,tp)  = lme2.Coefficients(:,8).Upper;
    brain(bvar).lowCI(:,tp) = lme2.Coefficients(:,7).Lower;   

end

%%% plot preSMA
subplot(1,2,bvar)
cIdx=1;
[~,ConIdx] = ismember('ic_1',lme2.CoefficientNames);
errBar = brain(bvar).upCI - brain(bvar).estim;
for i = ConIdx:ConIdx+3   
    hold on
    shadedErrorBar(time2plot,brain(bvar).estim(i,:),errBar(i,:),'lineProps',{'-','color',colorBeforIdx{cIdx}, 'LineWidth',1.5});
    xlim('tight'); 
    cIdx = cIdx +1;
end
legend(lme2.CoefficientNames(ConIdx:ConIdx+3))
ylim([-.25 .25])

%% plot average
figure(1); clf

ic_list = {'ACC' 'preSMA'};
typ_list = {'con' 'incon'};
before_list = {'incon' 'con'};
freq_list = {'congr50' 'congr20'};


bvar = 1;       
subplot(1,2,bvar)
list = [];
for typi=1:2   
    for condi=1:2

    trIdx = Tx.type==typ_list(typi) & Tx.freq==freq_list(condi); 
    sum(trIdx)
    hdrtyp = squeeze(mean(hdr2(trIdx,:,bvar))); 
    hdrerr = squeeze(std(hdr2(trIdx,:,bvar)))/sqrt(sum(trIdx)); 
    
    shadedErrorBar(time2plot,hdrtyp,hdrerr,'lineProps',{'-','color',colorFreqIdx{2*typi+condi-2}, 'LineWidth',1.5});
    hold on
    title(ic_list(bvar))
    list{end+1} = strcat(typ_list{typi}, '-', freq_list{condi}) ;
    end
end
legend(list)

bvar = 2;
subplot(1,2,bvar)
list = [];
for typi=1:2   
    for condi=1:2

    trIdx = Tx.type==typ_list(typi) & Tx.before==before_list(condi); 
    sum(trIdx)
    hdrtyp = squeeze(mean(hdr2(trIdx,:,bvar))); 
    hdrerr = squeeze(std(hdr2(trIdx,:,bvar)))/sqrt(sum(trIdx)); 
    
    shadedErrorBar(time2plot,hdrtyp,hdrerr,'lineProps',{'-','color',colorBeforIdx{2*typi+condi-2}, 'LineWidth',1.5});


    %shadedErrorBar(frex,psdrawsbAve(bl,:), psdrawsbSE(bl,:) ,'lineProps',{'-','color',colorIdx{bl}, 'LineWidth',1.5})

    
    hold on
    title(ic_list(bvar))
    list{end+1} = strcat(typ_list{typi}, '-p.', before_list{condi}) ;
    end
end
legend(list)
ylim([-.25 .25])

%%

bvar = 1;

tp = find(time2plot==4.5);

% extract brain value at tp
brainVal = nan(height(Tx),1);
for sb = sb_in
    brainVal(Tx.Subject==sb)  = squeeze(hdr2(Tx.Subject==sb,tp,bvar));    
end

% create new table with correct trials only
Tz = [Tx array2table(brainVal,'VariableNames',brainVars(bvar))]; 

%%% run LMM with behavioral variables
fitlme(Tz, strcat(brainVars{bvar}, '~  freq*type + normRT + (1 + freq*type + normRT|Subject)'))


%%
bvar = 2;
tp = find(time2plot==4);
% extract brain value at tp
brainVal = nan(height(Tx),1);
for sb = sb_in
    brainVal(Tx.Subject==sb)  = squeeze(hdr2(Tx.Subject==sb,tp,bvar));    
end

% create new table with correct trials only
Tz = [Tx array2table(brainVal,'VariableNames',brainVars(bvar))]; 

%%% run LMM with behavioral variables
fitlme(Tz, strcat(brainVars{bvar}, '~ before*type + normRT +  (1 + before*type + normRT |Subject)')) %
