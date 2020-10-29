%%% Temporal signatures of critical dynamics in cortical excitability as probed by early somatosensory responses %%%
% Authors: Tilman Stephani, Gunnar Waterstraat, Stefan Haufe, Garbiel Curio, Arno Villringer & Vadim V. Nikulin (2020)

%%% Code for DFA in thalamus-related CCA components %%%


%% prepare environment
clear

addpath('/data/p_01972/EEGLAB/eeglab14_1_1b/');
addpath('/data/p_01972/Scripts/Code_data_manuscript_N20_power_law_dynamics/Functions/')

eeglab
% set eeglab options to double precision
pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 1, 'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 0, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 0);


%% paths and names
savepath_pre = '/data/pt_01972/Preproc_data/N20_study1/'; % first part of save directory
datapath_pre = '/data/pt_01972/Preproc_data/N20_study1/'; % first part of data directory
savepath_figures = '/data/pt_01972/Results/figures/';
savepath_MANUSCRIPT = '/data/p_01972/Scripts/Code_data_manuscript_N20_power_law_dynamics/Data/'; % save path for preprocessed data that can be published according to the European GDPR

subj_dir = dir('/data/p_01972/RAW_EEG/N20_study1/*_*');
subj_names = {subj_dir.name};

% define conditions
cond = 'notask'; % passive stimulation sequence


%% Select subjects with good thalamus signal

% Subjects with thalamus-related CCA components and their polarities
% visually inspected in printed CCA plots (see N20_probe_manuscript_code_CCA_DFA.m)
subj_P15 = [ 2, 4, 5, 8, 9, 11, 12, 21, 22, 24, 25, 26, 29]; % subject
comp_P15 = [ 4, 3, 3, 3, 4,  3,  3,  3,  4,  4,  4,  3,  3]; % CCA component with thalamus activity
pol_P15  = [-1, 1, 1, 1, 1,  1,  1, -1, -1,  1,  1,  1,  1]; % polarity of P15 peak

load([savepath_pre 'Classify_CCA_notask_30to200Hz_nonotch_stdzd_CP4.mat']); % contains classification of other CCA components

% construct time vector
sr = 5000;
epoch_length_ms = [-100 600];
time_vec = linspace(epoch_length_ms(1), epoch_length_ms(2), 0.7*sr);


%% Calculate DFA time courses of thalamus-related CCA components
DFA_cca_thal = [];
DFA_cca_thal_perm = [];
mean_CCA_thalamus = [];
count = 0;
for s = subj_P15 % loop through all subjects that have a P15 component
    disp(subj_names{s})
    
    % load CCA data
    savepath_full = [savepath_pre subj_names{s} '/']; % complete save directory
    load([savepath_full 'CCA_' subj_names{s} '_' cond '_30to200Hz_nonotch_stdzd_CP4.mat'])

    count = count + 1; % subject counter    
   
    
    % DFA parameters
    start = 7; % smallest window size (trials) for fluctuation estimation
    stop = 70; % largest window size (trials) for fluctuation estimation
    num_segment=30; % number of windows
    start_fit = 7; % smallest window size (trials) for loglog regression line
    stop_fit  = 70; % largest window size (trials) for loglog regression line
    fg=0;  
    
    % prune parameters (account for excluded trials)
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
   
    % signals
    CCA_thal = squeeze(CCA_comps(comp_P15(count), :, :));
    CCA_thal_perm = CCA_thal(:,randperm(size(CCA_thal,2)));
    
    mean_CCA_thalamus(:,count) = mean(CCA_comps(comp_P15(count), :, :),3) * pol_P15(count); % average and standardize sign    
    
    % DFA loop
    parfor i = 1:size(CCA_comps, 2)
        [DFA_cca_thal(i,count),Amplitude,Alpha,time,st,epochs]=dfa_2018(CCA_thal(i,:), start,stop,num_segment,start_fit,stop_fit,prune, fg);
        [DFA_cca_thal_perm(i,count),Amplitude,Alpha,time,st,epochs]=dfa_2018(CCA_thal_perm(i,:),start,stop,num_segment,start_fit,stop_fit,prune, fg);
        
        if mod(i, 500) == 0, fprintf('sample #%d \n', i), end % show message every 100th trial
    end
    
end

% save DFA exponents
save([savepath_MANUSCRIPT 'DFA_exponents_time_thalamus_MANUSCRIPT.mat'], ...
    'DFA_cca_thal', 'DFA_cca_thal_perm', 'mean_CCA_thalamus', 'subj_P15', 'comp_P15', 'pol_P15')






