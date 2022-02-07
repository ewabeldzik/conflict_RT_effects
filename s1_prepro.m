% Script written 21.01.2020 by Ewa Beldzik
% EEG data preprocessing, part 1:
% (1) selecting time range and numbering trials
% (2) high & low pass filter 
% (3) average reference

clear; clc
pathdir = 'G:\SEF2\eeg\data\sef\';

for sb=1:37
% load datasets from BVA (all channels except for #32 - ECG)
EEG = pop_loadbv('G:\SEF2\eeg\raw\sefexp\',strcat('Expsef_', num2str(sb),'_1.vhdr'),[],[1:31,33:64]);

%%% (1) cutting resting periods; k-nb of trial within each block
trange=[]; k=1;
for i=1:size(EEG.event,2)
    EEG.event(i).trial=NaN;
    if strcmp(EEG.event(i).type,'S 99') 
    trange(j,1)= EEG.event(i).latency; 
    elseif strcmp(EEG.event(i).type,'S 88')
    trange(j,2)= EEG.event(i).latency;     
    j=j+1; k=1;
    elseif strcmp(EEG.event(i).type,'S 13')
    EEG.event(i).trial=k;
    k=k+1;
    end
end

if any(trange==0)
error(strcat('Missing info about time range for sb', num2str(sb)))
end

EEG = pop_select( EEG,'point', trange);
EEG = pop_selectevent( EEG, 'type',{'Sync On'},'select','inverse','deleteevents','on'); %'R' 'R128' 

% saving files00
EEG.setname=strcat('sb',num2str(sb),'-00');
if ~exist(strcat(pathdir,'sb',num2str(sb)),'dir'); mkdir(strcat(pathdir,'sb',num2str(sb)));end
EEG = pop_saveset( EEG, 'filename',strcat('sb',num2str(sb),'-00.set'),'filepath',strcat(pathdir,'sb',num2str(sb)));

% correct technical error (i.e. 80 ms screen delay)
for tr=1:length(EEG.event)
if strcmpi(EEG.event(tr).type,'S 13'); EEG.event(tr).latency=EEG.event(tr).latency+40; end
end

%%% (2) high & low pass filter 
EEG = pop_eegfiltnew(EEG, 'locutoff',0.5,'hicutoff',35,'plotfreqz',0);

%%% (3) average reference
EEG = pop_reref( EEG, []);

% saving files01
EEG.setname=strcat('sb',num2str(sb),'-01');
EEG = pop_saveset( EEG, 'filename',strcat('sb',num2str(sb),'-01.set'),'filepath',strcat(pathdir,'sb',num2str(sb)));

end