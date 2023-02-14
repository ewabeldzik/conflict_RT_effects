% Script written 25.07.2022 by Ewa Beldzik 
% conducts the trial selection analyses 
% aimed to verify the RT effect in EEG and fMRI 'conflict markers'
% 
% Results are presented in the Figure 3 in the manuscript entitled:
% "A thin line between conflict and reaction time effects on EEG and fMRI
% brain signals" (2023) by Ewa Beldzik & Markus Ullsperger


clear; clc

load('E:\SEF2\beh\beh_sef.mat');
load('E:\SEF2\eeg\matfiles\ic_mf_stim.mat')

% EEG time-frequncy analysis
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
        tf_temp(fi,:,:) = abs(as.^2) - mean(abs(as(bIdx,:).^2)); %./

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
timeBOLD = 0:0.1:9.9;
% define matrix for hdr & its peak
hdr = nan(height(T),length(timeBOLD),length(ic_list));
peak = nan(height(T),length(ic_list));
% time range for peak extraction
peakTime = timeBOLD > 3 & timeBOLD < 5;

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
            hdrSes(tr,:,:)=itc(tr_list(tr):tr_list(tr)+length(timeBOLD)-1,ic_list) - itc(tr_list(tr),ic_list);
        end
        
        hdr(SesIdx,:,:) = hdrSes;
        peak(SesIdx,:) = mean(hdrSes(:,peakTime,:),2);
         
    end

end


%% select trials based on RTs and plot
types = {'con','incon'};

% parameters for theta power extraction
thetaIdx = frex >= 4 & frex <= 8;
timeEEGIdx = time(2).t < 400; timeEEG = time(2).t(timeEEGIdx);
timeBresp = -300:0; timeBrespIdx = unique(dsearchn(time(2).t',timeBresp'));

% initialize output time-frequency data
brain = [];  
ic_brain  = labels(ic_list);


for sb = sb_in     
    
    brain(1).tc(sb,:,:) = squeeze(mean(tf(sb,2).data(thetaIdx,timeEEGIdx,:)));
    brain(1).peak(sb,:) = mean(brain(1).tc(sb,timeBrespIdx,:));
    
    brain(2).tc(sb,:,:) = squeeze(hdr(T.Subject==sb,:,find(ismember(ic_brain,'ACC'))))';
    brain(2).peak(sb,:) = squeeze(peak(T.Subject==sb,find(ismember(ic_brain,'ACC'))));
    
    brain(3).tc(sb,:,:) = squeeze(hdr(T.Subject==sb,:,find(ismember(ic_brain,'preSMA'))))';
    brain(3).peak(sb,:) = squeeze(peak(T.Subject==sb,find(ismember(ic_brain,'preSMA'))));

end    
   
%% plot the first column - 'traditional' congruency comparison
hFig = figure(1); clf
set( hFig, 'units','normalized','outerposition',[.1 .1 0.4 0.8]);

brain_list = {'theta', 'ACC','preSMA'};

% define colors
colorsStr = {'b','r'};
colors = [.34 .34 1; .9 .3 .3]; % blue & red
yLimits = [-.02 .08; -0.18 0.12; -0.1 .21]; 

data2plot = [];
rt=[]; how_many=[];
tval=[]; pval=[]; pcor = [];

% calculate mean
for sb = sb_in    

    % Calculating mean for each stimulus type
    Tsb = T(T.Subject==sb,:);

    for typ=1:length(types)     
            
        typIdx = Tsb.accur=='hit' & Tsb.badep==1 & Tsb.type==types{typ};        
        for bvar = 1:length(brain_list)
            rt(sb,typ) = mean(Tsb.stimRT(typIdx));
            how_many(sb,typ) = sum(typIdx);    
            brain(bvar).tctyp(typ,sb,:) = mean(brain(bvar).tc(sb,:,typIdx),3);               
        end

    end
end   
%
subplot(4,3,1)
for typ=1:length(types)
    h1 = raincloud_plot(rt(sb_in,typ), 'box_on', 1, 'color', colors(typ,:), 'alpha', 0.8,...
     'box_dodge', 1, 'box_dodge_amount', .2*typ, 'dot_dodge_amount', .2*typ,'box_col_match', 0);
    h1{1}.LineWidth = 1;
    h1{3}.LineWidth = 1;
    h1{5}.LineWidth = 1;
    h1{6}.LineWidth = 1;
    xlim([450 900]); xticks(500:100:800)
    hold on
end

for bvar = 1:length(brain_list)
    if bvar == 1;  time2plot = timeEEG; else; time2plot = timeBOLD; end
    
    brain2plot = [];
    [~,pval] = ttest(brain(bvar).tctyp(1,sb_in,:),brain(bvar).tctyp(2,sb_in,:));    
    pcor = 100*fdr_bh(squeeze(pval)); pcor(pcor==0)=NaN;  
    

    subplot(4,3,3*bvar+1)
    area(time2plot,pcor,-100,'LineStyle','none','FaceColor',[.87 .87 .87]);
    hold on
    for typ=1:length(types)
        brain2plot(typ,:) = nanmean(brain(bvar).tctyp(typ,sb_in,:));    
        errBar = nanstd(squeeze(brain(bvar).tctyp(typ,sb_in,:)))/sqrt(length(sb_in));     
        shadedErrorBar(time2plot,brain2plot(typ,:),errBar,'lineprops',colorsStr{typ});
    end
    if bvar==1;  xticks(-600:300:300); end
    ylim(yLimits(bvar,:)); xlim('tight')
end

[ mean(rt(sb_in,:)) mean(how_many(sb_in,:))]
[~, ps] = ttest(rt(sb_in,1),rt(sb_in,2))


% plot the second column - 'equalized RT' congruency comparison
data2plot = [];
rt=[]; how_many=[];
tval=[]; pval=[]; pcor = [];

% calculate mean
for sb = sb_in    

    % Calculating mean for each stimulus type
    Tsb = T(T.Subject==sb,:);

    for typ=1:length(types)     
            
        if strcmp(types{typ},'con')
            typIdx = Tsb.accur=='hit' & Tsb.badep==1 & Tsb.type==types{typ} & Tsb.normRT > -0.6 & Tsb.normRT < 1.3;   
        else
            typIdx = Tsb.accur=='hit' & Tsb.badep==1 & Tsb.type==types{typ} & Tsb.normRT > -1.5 & Tsb.normRT < 1.4;            
        end     
        rt(sb,typ) = mean(Tsb.stimRT(typIdx));
        how_many(sb,typ) = sum(typIdx);

        for bvar = 1:length(brain_list)    
            brain(bvar).tctyp(typ,sb,:) = mean(brain(bvar).tc(sb,:,typIdx),3);               
        end

    end
end   


subplot(4,3,2)
for typ=1:length(types)
    h1 = raincloud_plot(rt(sb_in,typ), 'box_on', 1, 'color', colors(typ,:), 'alpha', 0.8,...
     'box_dodge', 1, 'box_dodge_amount', .2*typ, 'dot_dodge_amount', .2*typ,'box_col_match', 0);
    hold on
    h1{1}.LineWidth = 1;
    h1{3}.LineWidth = 1;
    h1{5}.LineWidth = 1;
    h1{6}.LineWidth = 1;
    xlim([450 900]); xticks(500:100:800)
end

for bvar = 1:length(brain_list)
    if bvar == 1;  time2plot = timeEEG; else; time2plot = timeBOLD; end
    
    brain2plot = [];
    [~,pval] = ttest(brain(bvar).tctyp(1,sb_in,:),brain(bvar).tctyp(2,sb_in,:));    
    pcor = 100*fdr_bh(squeeze(pval)); pcor(pcor==0)=NaN;  
    
%   subplot(3,4,bvar+5)
    subplot(4,3,3*bvar+2)
    area(time2plot,pcor,-100,'LineStyle','none','FaceColor',[.87 .87 .87]);
    hold on
    for typ=1:length(types)
        brain2plot(typ,:) = nanmean(brain(bvar).tctyp(typ,sb_in,:));    
        errBar = nanstd(squeeze(brain(bvar).tctyp(typ,sb_in,:)))/sqrt(length(sb_in));     
        shadedErrorBar(time2plot ,brain2plot(typ,:),errBar,'lineprops',colorsStr{typ});
    end
    if bvar==1;  xticks(-600:300:300); end
    ylim(yLimits(bvar,:)); xlim('tight')
end

[ mean(rt(sb_in,:)) mean(how_many(sb_in,:))]
[~, ps] = ttest(rt(sb_in,1),rt(sb_in,2))



% plot the third row - 'equalized RT' for N bins 
binN = 4;

data2plot = []; rt2plot = []; tempBrain =[]; meanRT = [];
how_many = []; rt = [];
tval=[]; pval=[]; pcor = [];

corQuantiles = [.08 .04 .04 .08];

% normalize peak values for each sb
for sb=sb_in
    Tsb = T(T.Subject==sb,:);
    typIdx = (Tsb.type=='con' | Tsb.type=='incon') & Tsb.accur=='hit' & Tsb.badep==1 & Tsb.stimRT~=0;  

    for bvar = 1:3
        brain(bvar).peakZ(sb,typIdx) = zscore(brain(bvar).peak(sb,typIdx));
    end

end


for bvar = 1:length(brain_list)+1

    for typ = 1:length(types)
    for bin = 1:binN
        for sb=sb_in

            Tsb=T(T.Subject==sb,:);        
            Q = quantile(Tsb.normRT(Tsb.type=='incon'),0.05:0.9/binN:1); 
            
            if strcmp(types{typ},'con')
                typIdx = Tsb.accur=='hit' & Tsb.badep==1 & Tsb.type==types{typ} & Tsb.normRT > Q(bin) + corQuantiles(bin) & Tsb.normRT < Q(bin+1);   
            else
                typIdx = Tsb.accur=='hit' & Tsb.badep==1 & Tsb.type==types{typ} & Tsb.normRT > Q(bin) & Tsb.normRT < Q(bin+1);            
            end
            how_many(bin).data(sb,typ) = sum(typIdx);

            tempRT(sb) = mean(Tsb.stimRT(typIdx));

            if bvar > 1
                tempBrain(sb) = mean(brain(bvar-1).peakZ(sb,typIdx));
            end

        end

        if bvar == 1
        data2plot(bvar).data{bin,typ} = tempRT(sb_in);
        meanRT(typ,bin) = mean(tempRT(sb_in));
        else
        data2plot(bvar).data{bin,typ} = tempBrain(sb_in);
        end
    end
    end  


    % plot rainclouds
    subplot(4,3,3*bvar)

    h = rm_raincloud(data2plot(bvar).data, colors);

    for typ = 1:2
        for b = 1:binN
            h.s{b,typ}.SizeData = 5;
            h.m(b,typ).SizeData = 30;
            h.m(b,typ).MarkerEdgeAlpha = .8;
            h.m(b,typ).LineWidth = 1;
            if b < binN
                h.l(b,typ).Color = colors(typ,:)*0.6;
                h.l(b,typ).LineStyle = '--';
                h.l(b,typ).LineWidth = 1;
            end
        end
    end

    ylim('tight');  xlim('tight')
    
    for b = 1:binN
    [tval(bvar,b) pval(bvar,b)] = ttest(data2plot(bvar).data{b,1},data2plot(bvar).data{b,2});
    end

end
pval

