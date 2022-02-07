% Script written 10.10.2020 by Ewa Beldzik
% EEG data preprocessing, part 2
% Marking bad epochs

clear; clc
% load behavioural data (table T)
load('G:\SEF2\beh\beh_sef.mat'); 
pathdir = strcat('G:\SEF2\eeg\data\sef\');

N = 37;
ntrials = 1120;
badep = nan(height(T),1);

for sb=1:37
% load EEG datasets
EEG = pop_loadset('filename',strcat(pathdir,'sb',num2str(sb),'\sb',num2str(sb),'-01.set'));

%%% mark bad epochs based on variance during the stimulus/baseline period
EEGs = pop_epoch( EEG, {'S 13'}, [-0.5  1.2], 'epochinfo', 'yes');
% calculate mean variance for all channels
for i=1:EEGs.trials; epvar(i,1)=mean(var(EEGs.data(:,:,i),0,2)); end
epvar=zscore(epvar);

%%% mark bad epochs based on abs amplitude +/- 100ms around the response
EEGr = pop_epoch( EEG, {'S 30','S 32','S 46','S 48'}, [-1.2  0.6], 'epochinfo', 'yes');
% initialize output vector
respAmp = zeros(ntrials,1); T2=T(T.Subject==sb,:);  respIdx = T2.stimRT~=0;

% remove "double clicks", i.e. trials witout stimulus
if  sum(respIdx) ~= EEGr.trials

    for i=1:EEGr.trials
        rIdx = find(cell2mat(EEGr.epoch(i).eventlatency)==0);
        if ~ismember('S 13',EEGr.epoch(i).eventtype(1:rIdx))
            ep2rem = i;
            cprintf('red',['sb',num2str(sb),' - misclick in ',num2str(i),' trial\n'])
        else
            sIdx = find(ismember(EEGr.epoch(i).eventtype(1:rIdx),'S 13'));
            if cell2mat(EEGr.epoch(i).eventlatency(rIdx))-cell2mat(EEGr.epoch(i).eventlatency(sIdx))<20
                ep2rem = i;
                cprintf('red',['sb',num2str(sb),' - comision in ',num2str(i),' trial\n'])
            end
        end
    end
EEGr = pop_selectevent( EEGr, 'omitepoch',ep2rem);    
end

% calculate abs amp
EEGr = pop_epoch( EEGr, {'S 30','S 32','S 46','S 48'}, [-0.1  0.1], 'epochinfo', 'yes');
if  sum(respIdx) ~= EEGr.trials; cprintf('red',['sb',num2str(sb),' - check trials count']); end
respAmp(respIdx) = max(squeeze(max(EEGr.data,[],2)) - squeeze(min(EEGr.data,[],2)));

% index bad epochs
badep(T.Subject==sb,:)= (respAmp < 150) .* (abs(epvar) < 3);

% clear EEGs EEGr
end
T = [T table(badep)];
% save table manually (safer option;)