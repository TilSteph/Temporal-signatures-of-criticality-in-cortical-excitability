%%% Temporal signatures of critical dynamics in cortical excitability as probed by early somatosensory responses %%%
% Authors: Tilman Stephani, Gunnar Waterstraat, Stefan Haufe, Garbiel Curio, Arno Villringer & Vadim V. Nikulin.

%%% Code for obtaining pre-stimulus alpha envelopes from the same sources as the CCA components %%%



%% prepare environment
clear

addpath('/data/p_01972/EEGLAB/eeglab14_1_1b/');
addpath('/data/p_01972/Scripts/');
addpath('/data/p_01972/Scripts/Misc/export_fig-master');

eeglab

%% paths and names
savepath_pre = '/data/pt_01972/Preproc_data/N20_study1/'; % first part of save directory
datapath_pre = '/data/pt_01972/Preproc_data/N20_study1/'; % first part of data directory

subj_dir = dir('/data/p_01972/RAW_EEG/N20_study1/*_*');
subj_names = {subj_dir.name};

% define conditions
cond = 'notask';

%% subject loop
for s = 1:length(subj_names)
    if (strcmp(cond, 'task') && s==1) || (strcmp(cond, 'task') && s==13) % subject 13 has no task data; subject FC3 has been removed
        continue
    end
    
    subj = subj_names{s};
    fprintf(['\n\n ------------   SUBJECT ' subj '  condition: ' cond '   ------------ \n\n']);
    
    % add suffixes
    savepath_full = [savepath_pre subj '/']; % complete save directory
    data_path = [datapath_pre subj '/']; % ses file location
    subj = [subj '_' cond '_pchip'];
        
    % Load broad-band data  
    EEG = pop_loadset('filename',[subj '_sr5kHz_1to200Hz_vi_averef_nonotch_ICA_removed.set'],'filepath', [datapath_pre subj_names{s} '/']);
      
    % Epoch continuous data
    epoch_length_t = [-0.5 -0.005]; % in sec
    [EEG, ix_accepted_cont] = pop_epoch(EEG, {'A - Out', 'B - Out'}, epoch_length_t, 'newname', ['N20_probe_epochs'], 'epochinfo', 'yes'); % use second output later for choosing epochs (e.g. in behavioral data) AND prune parameters in DFA!!!!
    EEG.etc.accepted_epochs = ix_accepted_cont; % save for later
    
    % load previously calculated CCA data
    load([savepath_full 'CCA_' subj_names{s} '_' cond '_30to200Hz_nonotch_stdzd_CP4.mat'])    
    load([savepath_pre 'Classify_CCA_notask_30to200Hz_nonotch_stdzd_CP4.mat']) % corresponding classification of components
    
    
    %% trial loop
    filt_hz = [8 13]; % [4 7]; [8 10]; [10 13]; [15 20]; [20 30];
    [b,a] = butter(2, filt_hz/(EEG.srate/2));
    EEG_data_tmp =  zeros(2, EEG.pnts, EEG.trials); % initialize matrix
    for i = 1:size(EEG.data,3)        
        % band-pass filter prestimulus epochs
        sig = EEG.data(:,:,i);
        sig_mirr = [fliplr(sig) sig fliplr(sig)]; % mirror + padding (left side)
        sig_mirr = filtfilt(b, a, double(sig_mirr)'); % filter in alpha range 
        
        % apply CCA weights to prestimulus alpha
        sig_mirr_cca = Wst(:,[comps_tangential(s) comps_radial(s)])' * sig_mirr'; % only look at tangential and radial CCA components
        
        % get envelope of band-pass filtered data
        sig_mirr_cca_env = abs(hilbert(sig_mirr_cca'))'; 
        
        % take only middle part of mirrored signal
        EEG_data_tmp(:,:,i) = sig_mirr_cca_env(:, end/3+1:end/3*2); 
    end
    
    EEG.data = EEG_data_tmp; % put filtered signal (CCA components) back in EEG structure
    EEG = eeg_checkset(EEG); % correct other EEG fields    
    
    % save dataset
    EEG = pop_saveset(EEG, 'filename',[subj '_sr5kHz_' num2str(filt_hz(1)) 'to' num2str(filt_hz(2)) 'Hz_prestimulus_envelope_mirrored_nonotch_CCA.set'],'filepath', [datapath_pre subj_names{s} '/']);    
end
    

