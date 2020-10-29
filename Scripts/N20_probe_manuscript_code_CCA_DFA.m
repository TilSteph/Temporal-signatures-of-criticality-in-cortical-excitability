%%% Temporal signatures of critical dynamics in cortical excitability as probed by early somatosensory responses %%%
% Authors: Tilman Stephani, Gunnar Waterstraat, Stefan Haufe, Garbiel Curio, Arno Villringer & Vadim V. Nikulin.

%%% Code for single-trial analysis using CCA and DFA %%%


%% environment preparation
clear
addpath('/data/p_01972/EEGLAB/eeglab14_1_1b/');
addpath('/data/p_01972/Scripts/Code_data_manuscript_N20_power_law_dynamics/Functions/')
addpath('/data/p_01972/Scripts/Misc/export_fig-master'); % for exporting figures (from https://de.mathworks.com/matlabcentral/fileexchange/23629-export_fig)

eeglab

% set eeglab options to double precision (important for CCA)
pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 1, 'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 0, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 0);

% avoid "_" being interpreted as subscript
set(0, 'DefaulttextInterpreter', 'none')

%% paths and names
savepath_pre = '/data/pt_01972/Preproc_data/N20_study1/'; % first part of save directory
datapath_pre = '/data/pt_01972/Preproc_data/N20_study1/'; % first part of data directory
savepath_figures = '/data/pt_01972/Results/figures/';
savepath_MANUSCRIPT = '/data/p_01972/Scripts/Code_data_manuscript_N20_power_law_dynamics/Data/'; % save path for preprocessed data that can be published according to the European GDPR

subj_dir = dir('/data/p_01972/RAW_EEG/N20_study1/*_*');
subj_names = {subj_dir.name};

% define conditions
cond = 'notask'; % passive stimulation sequence


%% Subject loop
iCCA_tangential_s = [];
iCCA_radial_s = [];
for s = 1:33 %1:length(subj_names)   
    subj = subj_names{s};
    fprintf(['\n\n ------------   SUBJECT ' subj '  condition: ' cond '   ------------ \n\n']);
    
    % add suffixes
    savepath_full = [savepath_pre subj '/']; % complete save directory
    data_path = [datapath_pre subj '/']; % ses file location
    subj = [subj '_' cond '_pchip'];
    
    % Load data
    EEG = pop_loadset('filename',[subj '_sr5kHz_30to200Hz_vi_averef_nonotch_ICA_removed.set'],'filepath', [datapath_pre subj_names{s} '/']);
    
    % Convert data to double (according to eeglab options it should already be double)
    EEG.data = double(EEG.data);    
    
    % plot suffix (when saved)
    pl_suffx = ['_' cond '_30to200Hz_nonotch_stdzd_CP4']; % according to loaded data
    
    %% Prune parameters for later DFA (continuous segments without omitted trials)
    ix_break_start = [500 1200 1400 1600 1800];

    % adjust break indices to rejected segments
    EEG_marked = pop_loadset('filename',[subj_names{s} '_pchip_sr250Hz_1to45Hz_vi_notremoved.set'],'filepath', savepath_full); % this data also contains another experimental paradigm (see preprocessing script)

    ix_AB = find(ismember({EEG_marked.event.type}, {'A - Out', 'B - Out'})); 
    for j=1:length(ix_break_start), EEG_marked.event( ix_AB(ix_break_start(j))).type = 'last_stim'; end % rename last stimuli before breaks

    EEG_marked = pop_select(EEG_marked, 'point', [EEG_marked.etc.good_segments.start, EEG_marked.etc.good_segments.start+EEG_marked.etc.good_segments.dur]); % take only good segments

        % split datasets
        if strcmp(cond, 'notask')
            start_notask_tr = 11;
            end_notask_tr = 12;
            start_notask_lat = [EEG_marked.event(  find(ismember({EEG_marked.event.type}, num2str(start_notask_tr)), 1, 'last')  ).latency]; % get latency
            end_notask_lat = [EEG_marked.event(  find(ismember({EEG_marked.event.type}, num2str(end_notask_tr)), 1, 'last')  ).latency]; % get latency
            EEG_marked = pop_select(EEG_marked, 'point', [start_notask_lat, end_notask_lat]); % no task
        elseif strcmp(cond, 'task')        
            start_task_tr = 31;
            end_task_tr = 32;  
            start_task_lat = [EEG_marked.event(  find(ismember({EEG_marked.event.type}, num2str(start_task_tr)), 1, 'last')  ).latency]; % get latency
            end_task_lat = [EEG_marked.event(  find(ismember({EEG_marked.event.type}, num2str(end_task_tr)), 1, 'last')  ).latency]; % get latency
            EEG_marked = pop_select(EEG_marked, 'point', [start_task_lat, end_task_lat]); % task
        end

    % Count stimulation events before breaks
    ix_last_stim = find(ismember({EEG_marked.event.type}, {'last_stim'}));
    num_AB_before = []; % takes into account only breaks (not manually excluded data segments)
    for j = 1:length(ix_last_stim)
        num_AB_before(j) = length(find(ismember({EEG_marked.event(1:ix_last_stim(j)).type}, {'A - Out','B - Out'}))) + j; % find number of A-Out events before boundary; add j (last_stim events are stimulations, too)
    end

   

    
    %% Canonical Correlation Analysis (CCA)    
    % CCA parameters
    epoch_length_t = [-0.1 0.6]; % in sec
    epoch_markers = {'A - Out', 'B - Out'};
    train_win = [5 80]; % window used for training the CCA weights
    eval_win = epoch_length_t*1000;
    num_comp = 4;
    flag_standardize = 1; % or 0; if 1, first two components' polarity is standardized
    
    % Do CCA
    [CCA_comps, Wst, Ast, R_cca, iCCA_tangential, iCCA_radial, EEG, ix_accepted] = apply_CCA(EEG, epoch_markers, train_win, eval_win, num_comp, flag_standardize); % do CCA
    
    iCCA_tangential_s(s) = iCCA_tangential;
    iCCA_radial_s(s) = iCCA_radial;
    
    %% Adjustments needed because some epochs were too short
    % Take not-accepted epochs into account for DFA prune parameters
    ix_all_events = (1:length((ix_AB)))';
    ix_not_ex = ismember(ix_all_events, ix_accepted);
    ix_not_acc = find(ix_not_ex~=1);
    
    % correct previous num_AB_before
    for i = 1:length(num_AB_before)
        tmp_num_smaller = length(find(ix_not_acc < num_AB_before(i)));
        num_AB_before(i) = num_AB_before(i) - tmp_num_smaller;
    end
    
    
       
    %% save CCA-weighted data
    subj_cond = [subj_names{s} pl_suffx];
    if strcmp(cond, 'notask')
        save([savepath_full 'CCA_' subj_cond '.mat'], ...
            'CCA_comps', ...
            'Wst', ...
            'Ast', ... 
            'R_cca', ...
            'num_AB_before', ...
            'iCCA_tangential', ...
            'iCCA_radial', ...
            'subj_cond');
    elseif strcmp(cond, 'task')
        save([savepath_full 'CCA_' subj_cond '.mat'], ...
            'CCA_comps', ...
            'Wst', ...
            'Ast', ... 
            'R_cca', ...
            'num_AB_before', ...
            'iCCA_tangential', ...
            'iCCA_radial', ...
            'task_perf', ...
            'subj_cond');
    end 
    
    
    
    %% Plots of single subjects (these were visually inspected in order to classify the components as tangential, radial and thalamic components)
    title_str = {'CCA 1', 'CCA 2', 'CCA 3', 'CCA 4'};
    title_str{iCCA_tangential} = [title_str{iCCA_tangential} ' (tangential)'];
    title_str{iCCA_radial} = [title_str{iCCA_radial} ' (radial)'];
    
    % plot single trials
    num_comp_disp = 4; % show first four components
    x_lim = [-5 100];
    h = figure;
    set(h,'Visible','off');
    for k = 1:num_comp_disp
        subplot(4,1,k)
        imagesc(EEG.times, linspace(size(epochs_all, 2), 1, size(epochs_all, 2)), squeeze(CCA_comps(k,:,:))');
        colorbar; xlabel('time in ms'); ylabel('trials'); axis('xy')
        xlim(x_lim)
        caxis([-2 2])
        title(title_str{k})
    end
    suptitle_mod(['Subject ' subj_names{s} '  ' cond])
    
    export_fig([savepath_pre 'CCA_all_subjects' pl_suffx '.pdf'], '-pdf', '-append', h)
    
    
    % plot spatial patterns   
    ha = figure;
    set(ha, 'Visible', 'off');
    subplot(3,2,1:2); plot(R_cca,'.'), ylim([0 1])
    subplot(3,2,3); topoplot(Ast(:,1), EEG.chanlocs); title(title_str{1})
    subplot(3,2,4); topoplot(Ast(:,2), EEG.chanlocs); title(title_str{2})
    subplot(3,2,5); topoplot(Ast(:,3), EEG.chanlocs); title(title_str{3})
    subplot(3,2,6); topoplot(Ast(:,4), EEG.chanlocs); title(title_str{4})
    suptitle_mod(['CCA patterns  ' subj_names{s}])
    
    export_fig([savepath_pre 'CCA_patterns' pl_suffx '.pdf'], '-pdf', '-append', ha)
    
        
    % CCA average ERP
    num_comp_disp = 4; % show first four components
    x_lim = [-100 600];
    hh = figure;
    set(hh, 'Visible', 'off')
    for k = 1:num_comp_disp
        subplot(2,2,k)
        plot(EEG.times, mean(CCA_comps(k,:,:), 3))
        xlabel('time (ms)')
        ylabel('component activation (a.u.)')
        xlim(x_lim)
        title(title_str{k})
    end
    suptitle_mod(['CCA averages  ' subj_names{s}])

    export_fig([savepath_pre 'CCA_average' pl_suffx '.pdf'], '-pdf', '-append', hh)
     
end

    
    
%% Classify CCA components (according to visual inspection of time courses and spatial patterns)

comps_tangential = [2, 3, 2, 1, 1, 2, 2, 1, 1, 2, ...
                    1, 1, 1, 2, 1, 1, 1, 2, 2, 2, ...
                    1, 1, 2, 1, 2, 1, 2, 2, 1, 2, ...
                    1, 1, 1];

comps_radial     = [1, 1, 1, 2, 2, 1, 1, 2, 2, 1, ...
                    2, 2, 1, 1, 2, 2, 2, 1, 1, 1, ...
                    2, 1, 1, 2, 1, 2, 1, 1, 2, 1, ...
                    2, 2, 2];
save([savepath_pre 'Classify_CCA_notask_30to200Hz_nonotch_stdzd_CP4.mat'], 'comps_tangential', 'comps_radial')
       
                


%% Detrended Fluctuation Analysis (DFA) + save mean of CCA over all participants
cond = 'notask';
load([savepath_pre 'Classify_CCA_notask_30to200Hz_nonotch_stdzd_CP4.mat']) 

CCA_average = [];
CCA_tangential_average = [];
CCA_radial_average = [];
Ast_all = [];
Ast_tangential = [];
Ast_radial = [];

DFA_tangential = zeros(size(CCA_comps,2), length(subj_names));
DFA_radial = zeros(size(CCA_comps,2), length(subj_names));
R2_tang = zeros(size(CCA_comps,2), length(subj_names)); R2_rad = zeros(size(CCA_comps,2), length(subj_names));

n_shuffle_iter = 1000; % iterations of shuffling the trial order
DFA_tangential_perm = zeros(size(CCA_comps,2), length(subj_names), n_shuffle_iter);

for s = 1:length(subj_names)   
    subj = subj_names{s};
    fprintf(['\n\n ------------ Subject ' subj '  condition: ' cond '   ------------ \n\n']);
   
    % load CCA data  
    pl_suffx = ['_' cond '_30to200Hz_nonotch_stdzd_CP4']; % according to loaded data
    savepath_full = [savepath_pre subj '/']; % complete save directory
    load([savepath_full 'CCA_' subj_names{s} pl_suffx '.mat'])

    %% Average CCA components
    for k = 1:4 % loop through first four CCA components
        CCA_average(k, :, s) = mean(CCA_comps(k, :, :),3);
    end
    
    CCA_tangential_average(:, s) = mean(CCA_comps(comps_tangential(s), :, :),3);
    CCA_radial_average(:, s) = mean(CCA_comps(comps_radial(s), :, :),3);
    
    
    %% Spatial patterns
    for k = 1:4 % loop through first four CCA components
        Ast_all(k,:,s) = Ast(:,k); % 01_BT needs to have the same channel number for this!
    end
    
    Ast_tangential(:,s) = Ast(:, comps_tangential(s));
    Ast_radial(:,s) = Ast(:, comps_radial(s));
     
    %% DFAs over time    
    % DFA parameters
    start = 7; % smallest window size (trials) for fluctuation estimation
    stop = 70; % largest window size (trials) for fluctuation estimation
    num_segment=30; % number of windows
    start_fit = 7; % smallest window size (trials) for loglog regression line
    stop_fit  = 70; % largest window size (trials) for loglog regression line
    fg=0; 
    
    % prune parameters (do not calculate DFA across breaks in the data (i.e., excluded trials))
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
    
    % Do DFA    
    parfor i = 1:size(CCA_comps, 2)
        [DFA_tangential(i,s),Amplitude,Alpha,time,R2_tang(i,s),epochs]=dfa_2018(squeeze(CCA_comps(comps_tangential(s), i, :)),start,stop,num_segment,start_fit,stop_fit,prune, fg);
        [DFA_radial(i,s),Amplitude,Alpha,time,R2_rad(i,s),epochs]=dfa_2018(squeeze(CCA_comps(comps_radial(s), i, :)),start,stop,num_segment,start_fit,stop_fit,prune, fg);        
            
        % permute trials and calculate DFA
        for k = 1:n_shuffle_iter
            CCA_tang_perm = CCA_comps(comps_tangential(s), i, randperm(size(CCA_comps,3))); % permute trial order ("null hypothesis")
            [DFA_tangential_perm(i,s,k),Amplitude,Alpha,time,st,epochs]=dfa_2018(squeeze(CCA_tang_perm),start,stop,num_segment,start_fit,stop_fit,prune, fg);            
        end
        
        if mod(i, 500) == 0, fprintf('sample #%d \n', i), end % show message every 500th trial
    end
    
%     % Do CCA including variances for all window sizes (in order to be able to evaluate model fit)
%     Alpha = zeros(size(CCA_comps, 2),num_segment);
%     parfor i = 1:size(CCA_comps, 2)
%         [DFA_tangential(i,s),Amplitude,Alpha(i,:),time,st,epochs]=dfa_2018(squeeze(CCA_comps(comps_tangential(s), i, :)),start,stop,num_segment,start_fit,stop_fit,prune, fg);
%         %[DFA_radial(i,s),Amplitude,Alpha,time,st,epochs]=dfa_2018(squeeze(CCA_comps(comps_radial(s), i, :)),start,stop,num_segment,start_fit,stop_fit,prune, fg);                
%         if mod(i, 500) == 0, fprintf('sample #%d \n', i), end % show message every 100th trial
%     end
%     save([savepath_full subj '_notask_DFA_7to70_variance_per_window.mat'],'Alpha', 'time')
end

% save grand-averages and DFA results
save([savepath_MANUSCRIPT 'DFA_exponents_with_permutations_and_CCA_averages_MANUSCRIPT.mat'], ...
    'CCA_average', 'CCA_tangential_average', 'CCA_radial_average', 'Ast_all', 'Ast_tangential', 'Ast_radial', ...
    'DFA_tangential', 'DFA_radial', 'DFA_tangential_perm', 'R2_tang', 'R2_rad');






