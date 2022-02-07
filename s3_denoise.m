% Script written 10.10.2020 by Ewa Beldzik 
% EEG data preprocessing, part 3
% ICA denoising

clear; clc
% load behavioural data (table T)
load('G:\SEF2\beh\beh_sef.mat'); 
pathdir = strcat('G:\SEF2\eeg\data\sef\');

for sb=1:37
% load EEG datasets
EEG = pop_loadset('filename',strcat(pathdir,'sb',num2str(sb),'\sb',num2str(sb),'-01.set'));

% mark bad epochs


% high & low pass filter 
EEGtemp = pop_eegfiltnew(EEG, 'locutoff',4,'hicutoff',8,'plotfreqz',0);

% create temporary file - extract epochs 
EEGtemp = pop_epoch( EEGtemp, {'S 13'}, [0  1.2], 'epochinfo', 'yes');   

% select epochs of temporary file
T2=T(T.Subject==sb,:); 
epochIdx =(T2.type=='con' | T2.type=='incon') & T2.stimACC==1 & T2.badep==1;    
EEGtemp = pop_select(EEGtemp, 'trial', find(epochIdx==1));

% run ICA on temporary file
EEGtemp = pop_runica(EEGtemp, 'icatype', 'runica', 'extended',1,'interrupt','off');

% overwrite IC decomposition
EEG.icasphere  = EEGtemp.icasphere;
EEG.icaweights = EEGtemp.icaweights;
EEG.icawinv    = EEGtemp.icawinv;
EEG.icachansind= EEGtemp.icachansind;

% saving files01
EEG = pop_saveset( EEG, 'filename',strcat('sb',num2str(sb),'-01.set'),'filepath',strcat(pathdir,'sb',num2str(sb)));   

% epoching and saving files02
EEG = pop_epoch( EEG, {'S 13'}, [-1  1.8], 'epochinfo', 'yes');
EEG.setname=strcat('sb',num2str(sb),'-02-ica');
EEG = pop_saveset( EEG, 'filename',strcat('sb',num2str(sb),'-02-ica.set'),'filepath',strcat(pathdir,'sb',num2str(sb)));   

end

%% Initialiaze matrieces with frontocentral data
clear; clc
ic_nb = []; % which IC
ic_tc = []; % time-series of selected IC
ic_topo = []; % topography of selected IC

%% Next, plot components for each subject 

% sb = 1;
% EEG = pop_loadset('filename',strcat(pathdir,'sb',num2str(sb),'\sb',num2str(sb),'-02-ica.set'));
% EEG.icaact = (EEG.icaweights*EEG.icasphere)*reshape(EEG.data(EEG.icachansind,:), EEG.nbchan, EEG.trials*EEG.pnts);
% EEG.icaact = reshape(EEG.icaact, size(EEG.icaact,1),EEG.pnts, EEG.trials);
% pop_topoplot(EEG, 0, 1:63 ,strcat('sb',num2str(sb)),[] ,0,'electrodes','on');

%% Mark the frontocentral IC

% ic_nb(sb) = ;

%% Save IC time-course - for stimulus-locked analyses

% ic_tc(sb).data = squeeze(EEG.icaact(ic_nb(sb),:,:));
% ic_topo(sb,:) = EEG.icawinv(:,ic_nb(sb))*sign(EEG.icawinv(17,ic_nb(sb)));

%% After all sb, save data manually (safer option;) with one EEG file (for parameters)

% EEG.data = [];
% save as 'ic_mf_stim.mat'

%% Save IC time-course - for continuous analyses

% sb_in = find(ic_nb~=0); 
% 
% %%% initialize new matrix for IC time-series & session onset info
% ic_tc = []; ses_onset = [];
% 
% for sb = sb_in
%     EEG = pop_loadset('filename',strcat(pathdir,'sb',num2str(sb),'\sb',num2str(sb),'-01.set'));
%     EEG.icaact = (EEG.icaweights*EEG.icasphere)*reshape(EEG.data(EEG.icachansind,:), EEG.nbchan, EEG.trials*EEG.pnts);
%     ic_tc(sb).data = EEG.icaact(ic_nb(sb),:);
% 
%     ses=1;
%     for i=1:size(EEG.event,2)
%         if strcmp(EEG.event(i).type,'S 99')
%            ses_onset(sb,ses)=EEG.event(i).latency;
%            ses=ses+1;
%         end        
%     end
%     ses_onset(sb,6)=EEG.pnts;
% end

%% After all sb, save data manually (safer option;) with one EEG file (for parameters)

% save as 'ic_mf_contin.mat'