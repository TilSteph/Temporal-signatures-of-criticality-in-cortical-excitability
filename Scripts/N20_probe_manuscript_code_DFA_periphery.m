%%% Temporal signatures of critical dynamics in cortical excitability as probed by early somatosensory responses %%%
% Authors: Tilman Stephani, Gunnar Waterstraat, Stefan Haufe, Garbiel Curio, Arno Villringer & Vadim V. Nikulin.

%%% Code for DFA in CNAP of the median nerve (peripheral data) %%%


%% prepare environment
clear

addpath('/data/p_01972/EEGLAB/eeglab14_1_1b/');
addpath('/data/p_01972/Scripts/Code_data_manuscript_N20_power_law_dynamics/Functions/')

eeglab


%% paths and names
savepath_pre = '/data/pt_01972/Preproc_data/N20_study1/'; % first part of save directory
datapath_pre = '/data/pt_01972/Preproc_data/N20_study1/'; % first part of data directory
savepath_MANUSCRIPT = '/data/p_01972/Scripts/Code_data_manuscript_N20_power_law_dynamics/Data/'; % save path for preprocessed data that can be published according to the European GDPR

subj_dir = dir('/data/p_01972/RAW_EEG/N20_study1/*_*');
subj_names = {subj_dir.name};

% define conditions
cond = 'notask'; % passive stimulation sequence
    
DFA_exponents_time_peri = [];
mean_SEP = [];
%% subject loop
for s = 1:length(subj_names)

    subj = subj_names{s};
    fprintf(['\n\n ------------ generate and save DFA time courses of subject ' subj '  condition: periphery ' cond '   ------------ \n\n']);

    % load CCA data    
    savepath_full = [savepath_pre subj '/']; % complete save directory
    load([savepath_full 'CCA_' subj_names{s} '_' cond '_30to200Hz_stdzd.mat']) % only to obtain prune parameters fro DFA

    EEG = pop_loadset('filename',[subj_names{s} '_peri_notask_pchip_sr5kHz_HP70Hz_vi_notch.set'],'filepath', savepath_full);
    
    % Cut data into epochs
    epoch_length_t = [-0.1 0.6]; % in sec; shoudl be same length as where prune parameters are from
    [EEG, ix_accepted_events] = pop_epoch(EEG, {'A - Out', 'B - Out'}, epoch_length_t, 'newname', ['N20_probe_epochs'], 'epochinfo', 'yes'); % use second output later for choosing epochs (e.g. in behavioral data) AND prune parameters in DFA!!!!
        
    if s ~= 29
        mean_SEP(:,s) = mean(EEG.data(1,:,:),3); % average CNAP (median nerve)
    elseif s == 29
        mean_SEP(:,s) = mean(EEG.data(2,:,:),3); % second channel contains CNAP for this subject
    end
    
    % prune parameters
    n_epochs_all = size(CCA_comps, 3);
    prune_ends = [num_AB_before, n_epochs_all]; % prune_ends = [num_A_before, length(p1_epochs_A_B)];
    prune_len = []; % lengths of uninterrupted segments
    for i = 1:length(prune_ends)
        if i==1
            prune_len(i) = prune_ends(i)-1;
        else
            prune_len(i) = prune_ends(i)-prune_ends(i-1)-1;
        end
    end
    prune_starts = prune_ends - prune_len;
    prune = [prune_starts; prune_len]';

    % DFA parameters
    start = 7; % smallest window size (trials) for fluctuation estimation
    stop = 70; % largest window size (trials) for fluctuation estimation
    num_segment=30; % number of windows
    start_fit = 7; % smallest window size (trials) for loglog regression line
    stop_fit  = 70; % largest window size (trials) for loglog regression line
    fg=0;  

    % DFA exponents for time course
    if s ~= 29
        sig1 = squeeze(EEG.data(1,:,:));
    elseif s == 29
        sig1 = squeeze(EEG.data(2,:,:));
    end
    
    exponent1 = [];       
    parfor i = 1:size(EEG.data, 2)          
        [exponent1(i),Amplitude,Alpha,time,st,epochs]=dfa_2018(sig1(i,:),start,stop,num_segment,start_fit,stop_fit,prune, fg);
        if mod(i, 500) == 0, fprintf('sample #%d \n', i), end % show message every 500th sample
    end
    
    DFA_exponents_time_peri(:, s) = exponent1;
end

save([savepath_MANUSCRIPT 'DFA_exponents_time_CNAP_MANUSCRIPT.mat'], 'DFA_exponents_time_peri', 'mean_SEP')







 