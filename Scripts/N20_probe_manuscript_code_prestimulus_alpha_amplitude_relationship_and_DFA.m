%%% Temporal signatures of critical dynamics in cortical excitability as probed by early somatosensory responses %%%
% Authors: Tilman Stephani, Gunnar Waterstraat, Stefan Haufe, Garbiel Curio, Arno Villringer & Vadim V. Nikulin (2020)

%%% Code for examination of relationship between pre-stimulus alpha activity and SEP both regarding their amplitudes and DFA exponents %%%


%% prepare environment
clear

addpath('/data/p_01972/EEGLAB/eeglab14_1_1b/');
addpath('/data/p_01972/Scripts/Code_data_manuscript_N20_power_law_dynamics/Functions/')

eeglab


%% paths and names
savepath_pre = '/data/pt_01972/Preproc_data/N20_study1/'; % first part of save directory
datapath_pre = '/data/pt_01972/Preproc_data/N20_study1/'; % first part of data directory

subj_dir = dir('/data/p_01972/RAW_EEG/N20_study1/*_*');
subj_names = {subj_dir.name};

% define conditions
cond = 'notask';


%% Amplitude relationship
group_var = struct();
nbins = 5; % binning approach for visualization
TangCCA_bins_mean = [];

for s = 1:33
    disp(subj_names{s})
    savepath_full = [savepath_pre subj_names{s} '/']; % complete save directory
    subj = [subj_names{s} '_' cond '_pchip'];

    % load prestimulus alpha data
    filt_hz = [8 13];
    EEG = pop_loadset('filename',[subj '_sr5kHz_' num2str(filt_hz(1)) 'to' num2str(filt_hz(2)) 'Hz_prestimulus_envelope_mirrored_nonotch_CCA.set'],'filepath', savepath_full); % obtained from Get_prestimulus_alpha_cca_filtered.m

    % load previously calculated CCA data
    load([savepath_full 'CCA_' subj_names{s} '_' cond '_30to200Hz_nonotch_stdzd_CP4.mat'])    
    load([savepath_pre 'Classify_CCA_notask_30to200Hz_nonotch_stdzd_CP4.mat']) % classification of components

    % load continuous data
    EEG4CCA = pop_loadset('filename', [subj_names{s} '_' cond '_pchip_sr5kHz_30to200Hz_vi_averef_nonotch_ICA_removed.set'],'filepath', savepath_full);
    eval_win = [-100 600];
    [EEG4CCA, ix_accepted_cca] = pop_epoch(EEG4CCA, {'A - Out', 'B - Out'}, eval_win/1000, 'epochinfo', 'yes'); % use second output later for choosing epochs (e.g. in behavioral data) AND prune parameters in DFA
       
    % take only overlapping epochs
    [ix_overlap_alpha] = ismember(EEG.etc.accepted_epochs, ix_accepted_cca);
    [ix_overlap_cca] = ismember(ix_accepted_cca, EEG.etc.accepted_epochs);

    EEG.data = EEG.data(:,:,ix_overlap_alpha);
    CCA_comps = CCA_comps(:,:,ix_overlap_cca);

    % average prestimulus alpha in time bins
    epoch_length_t = [-0.5 -0.005];
    win_prestim_ms = [-200 -10]; % pre-stimulus extraction window          
    win_prestim_pt = (win_prestim_ms - epoch_length_t(1)*1000) .* EEG.srate/1000;

    alpha_prestim_cca = squeeze(mean(EEG.data(:,win_prestim_pt(1):win_prestim_pt(2),:),2)); % average prestim alpha; all CCA components

    % find out individual N20 peak
    epoch_length_SEP_t = [-0.1 0.6];
    mean_CCA_tangential = mean(CCA_comps(comps_tangential(s),:,:),3);
    N20_win = ([17 25] - epoch_length_SEP_t(1)*1000) .* EEG.srate/1000; % raw window for N20 search
    [~, peak_lat_raw] = min(mean_CCA_tangential(N20_win(1):N20_win(2)));
    peak_lat = peak_lat_raw + N20_win(1);
    
    % extract amplitude information from SEP        
    N20_win = [peak_lat-10 peak_lat+10]; % subject-specific window 2 ms around average peak

    % find peak amplitude
    [peak_amp_N20, peak_lat_raw] = min(squeeze(CCA_comps(comps_tangential(s), N20_win(1):N20_win(2), :)));
    if s == 21, [peak_amp_N20, peak_lat_raw] = min(squeeze(-1 * CCA_comps(comps_tangential(s), N20_win(1):N20_win(2), :))); end % flip polarity in subject 21 (automatic standardization did not work)
    
    % put in group variable
    group_var(s).alpha_prestim_cca = alpha_prestim_cca; % corresponding to tangential and radial CCA components
    group_var(s).peak_amp_N20 = peak_amp_N20;        
   
    % get SEP (tangential CCA) binned by prestim alpha (for visualization)
    [~, sort_ix] = sort(alpha_prestim_cca(1, :), 'ascend'); % tangential CCA
    ntrials = length(alpha_prestim_cca);
    for n = 1:nbins
        TangCCA_bins_mean(:,n,s) = mean(CCA_comps(comps_tangential(s), :, sort_ix(round(ntrials/nbins*(n-1))+1:round(ntrials/nbins*n)) ),3);
    end
end

% save
save([savepath_pre 'Predict_N20_from_prestim_alpha_200to10ms_CCA_filter_MANUSCRIPT.mat'], 'group_var', 'TangCCA_bins_mean')


% calculate correlation 
r_alpha_SEP = []; p_alpha_SEP = []; r_alpha_SEP_sensor = []; p_alpha_SEP_sensor = []; r_alpha_Mean_amp_SEP = []; p_alpha_Mean_amp_SEP = []; r_alpha_peak_amp_SEP = []; p_alpha_peak_amp_SEP = [];
alpha_prestim_cca_concat = []; peak_amp_N20_concat = []; r_alpha_sensor_peak_amp = []; p_alpha_sensor_peak_amp = [];
corr_type = 'Spearman';
for s = 1:33
    [r_alpha_peak_amp_SEP(s), p_alpha_peak_amp_SEP(s)] = corr(group_var(s).alpha_prestim_cca(1,:)', group_var(s).peak_amp_N20', 'type', corr_type); % tangential CCA component    

    alpha_prestim_cca_concat = [alpha_prestim_cca_concat, group_var(s).alpha_prestim_cca(1,:)];
    peak_amp_N20_concat = [peak_amp_N20_concat, group_var(s).peak_amp_N20];
end



% t null distribution on group level using surrogates
subj_in = [1:12,14:17,19:33];
r_alpha_peak_amp_SEP_z = atanh(r_alpha_peak_amp_SEP); % Fisher z transform
[H,P,CI,STATS] = ttest(r_alpha_peak_amp_SEP_z(subj_in)); % empirical data

t_emp = STATS.tstat; % empirical t value
nreps = 10000;
r_surr_z_mean = zeros(nreps, 1);
p_param_surr = zeros(nreps, 1);
for k = 1:nreps    
    subj_count = 0;
    for s = subj_in % loop through subjects
        subj_count = subj_count + 1;
        prestim_surrogate = AAFT_surrogate(group_var(s).alpha_prestim_cca(1,:), 1); % generate one surrogate time series
        r_tmp(subj_count) = corr(prestim_surrogate, group_var(s).peak_amp_N20', 'type', 'Spearman');
    end
    r_surr_z = atanh(r_tmp);
    r_surr_z_mean(k) = tanh(mean(r_surr_z)); % surrogate group-level correlation   
    
    [~,p_param_surr(k),~,STATS_surr] = ttest(r_surr_z); % surrogate data
    %[~,p_param_surr(k),~,STATS_surr] = ttest(r_surr_z, r_surr_z_mean); % compare one surrogate dataset against the average null correlation
    t_surr(k) = STATS_surr.tstat; % surrogate t value
    
    if mod(k, 1000)==0, display(['permutation #' num2str(k)]), end
end

p_val = mean(abs(r_surr_z_mean) >= abs(mean(r_alpha_peak_amp_SEP_z))); % p value from permutation test
mean_r = tanh(mean(r_alpha_peak_amp_SEP_z)); % average correlation



% prepare data for LME models to be computed in R
subj_in = [1:12,14:17,19:33];

alpha_prestim_cca_concat = []; peak_amp_N20_concat = []; subj_count = 0; subj_id = [];
for s = subj_in
    subj_count = subj_count + 1;
    alpha_prestim_cca_concat = [alpha_prestim_cca_concat, group_var(s).alpha_prestim_cca];
    peak_amp_N20_concat = [peak_amp_N20_concat, group_var(s).peak_amp_N20];
    subj_id = [subj_id; repmat(subj_count, length(group_var(s).peak_amp_N20), 1)];
end

alpha_prestim_tang_concat = alpha_prestim_cca_concat(1,:); % only tangential CCA component

save([savepath_pre 'LME_prestim_peak_amp_MANUSCRIPT'], 'alpha_prestim_tang_concat', 'peak_amp_N20_concat', 'subj_id');



%% DFA exponents of pre-stimulus alpha activity
subj_dir = dir('/data/p_01972/RAW_EEG/N20_study1/*_*');
subj_names = {subj_dir.name};

% define conditions
cond = 'notask';

nreps = 1000; % repetitions for DFA surrogates
DFA_alpha_prestim_tang = []; DFA_alpha_prestim_rad = []; DFA_alpha_prestim_tang_perm = zeros(nreps,length(subj_names));
for s = 1:length(subj_names)
    disp(subj_names{s})
    savepath_full = [savepath_pre subj_names{s} '/']; % complete save directory
    subj = [subj_names{s} '_' cond '_pchip'];

    % load prestimulus alpha data
    filt_hz = [8 13];
    EEG = pop_loadset('filename',[subj '_sr5kHz_' num2str(filt_hz(1)) 'to' num2str(filt_hz(2)) 'Hz_prestimulus_envelope_mirrored_nonotch_CCA.set'],'filepath', savepath_full); % cut until -5 ms

    % load previously calculated CCA data
    load([savepath_full 'CCA_' subj_names{s} '_' cond '_30to200Hz_nonotch_stdzd_CP4.mat'])    
    load([savepath_pre 'Classify_CCA_notask_30to200Hz_nonotch_stdzd_CP4.mat']) % classification of components

    % load continuous data
    EEG4CCA = pop_loadset('filename', [subj_names{s} '_' cond '_pchip_sr5kHz_30to200Hz_vi_averef_nonotch_ICA_removed.set'],'filepath', savepath_full);
    eval_win = [-100 600];
    [EEG4CCA, ix_accepted_cca] = pop_epoch(EEG4CCA, {'A - Out', 'B - Out'}, eval_win/1000, 'epochinfo', 'yes'); % use second output later for choosing epochs (e.g. in behavioral data) AND prune parameters in DFA!!!!
       
    % take only overlapping epochs
    [ix_overlap_alpha] = ismember(EEG.etc.accepted_epochs, ix_accepted_cca);
    [ix_overlap_cca] = ismember(ix_accepted_cca, EEG.etc.accepted_epochs);

    EEG.data = EEG.data(:,:,ix_overlap_alpha);
    CCA_comps = CCA_comps(:,:,ix_overlap_cca);

    % average prestimulus alpha in time bins
    epoch_length_t = [-0.5 -0.005];
    win_prestim_ms = [-200 -10]; % pre-stimulus extraction window           XXXXXX set parameter!
    win_prestim_pt = (win_prestim_ms - epoch_length_t(1)*1000) .* EEG.srate/1000;

    alpha_prestim_cca = squeeze(mean(EEG.data(:,win_prestim_pt(1):win_prestim_pt(2),:),2)); % average prestim alpha; all CCA components        
    
    % DFA on mean pre-stimulus power
    % parameters
    start = 7; % smallest window size (trials) for fluctuation estimation
    stop = 70; % largest window size (trials) for fluctuation estimation
    num_segment=30; % number of windows
    start_fit = 7; % smallest window size (trials) for loglog regression line
    stop_fit  = 70; % largest window size (trials) for loglog regression line
    fg=0;  
    
    % prune (account for excluded trials)
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
    
    
    % tangential
    [DFA_alpha_prestim_tang(s),Amplitude,Alpha,time,st,epochs]=dfa_2018(alpha_prestim_cca(1,:),start,stop,num_segment,start_fit,stop_fit,prune, fg); % with prune parameters
    
    % radial
    [DFA_alpha_prestim_rad(s),Amplitude,Alpha,time,st,epochs]=dfa_2018(alpha_prestim_cca(2,:),start,stop,num_segment,start_fit,stop_fit,prune, fg); % with prune parameters
        
    % permuted tangential CCA
    parfor k = 1:nreps
        alpha_prestim_cca_perm = alpha_prestim_cca(:,randperm(size(alpha_prestim_cca,2))); % shuffle trial order
        [DFA_alpha_prestim_tang_perm(k,s), Amplitude,Alpha,time,st,epochs]=dfa_2018(alpha_prestim_cca_perm(1,:),start,stop,num_segment,start_fit,stop_fit,prune, fg); % with prune parameters
    end
end

% save
save([savepath_pre 'DFA_timecourse_prestim_alpha_cca_MANUSCRIPT.mat'], 'DFA_alpha_prestim_tang', 'DFA_alpha_prestim_rad', 'DFA_alpha_prestim_tang_perm')


% statistical test against random DFA
subj_in = [1:12,14:17,19:33];
mean(DFA_alpha_prestim_tang(subj_in))
mean(mean(DFA_alpha_prestim_tang_perm(:,subj_in),1))

[H, p, CI, Stats] = ttest(DFA_alpha_prestim_tang(subj_in), mean(mean(DFA_alpha_prestim_tang_perm(:,subj_in),1))); % significant


%% DFA relationship between prestimulus alpha and SEP
load([savepath_pre 'DFA_timecourse_prestim_alpha_cca_MANUSCRIPT.mat']) % DFA exponents of pre-stimulus alpha activity
load([savepath_pre 'DFA_exponents_with_permutations_and_CCA_averages_MANUSCRIPT.mat']); % DFA exponent time course of CCA components

subj_in = [1:12, 14:17, 19:33]; % select subjects
srate = 5000; % sampling rate

% extract rms of DFA exponent time course
epoch_length_t = [-0.1 0.6];
win1_signal_ms = [20 25] - epoch_length_t(1)*1000; % in ms
win1_signal = win1_signal_ms * srate/1000;

win2_signal_ms = [25 30] - epoch_length_t(1)*1000; % in ms
win2_signal = win2_signal_ms * srate/1000;

win3_signal_ms = [30 35] - epoch_length_t(1)*1000; % in ms
win3_signal = win3_signal_ms * srate/1000;

win4_signal_ms = [35 40] - epoch_length_t(1)*1000; % in ms
win4_signal = win4_signal_ms * srate/1000;

mean1_DFA = []; mean2_DFA = []; mean3_DFA = []; mean4_DFA = [];
for s = 1:33    
    % arithmetic mean in time window of interest
    mean1_DFA(s) = mean(DFA_tangential(win1_signal(1):win1_signal(2), s));
    mean2_DFA(s) = mean(DFA_tangential(win2_signal(1):win2_signal(2), s));   
    mean3_DFA(s) = mean(DFA_tangential(win3_signal(1):win3_signal(2), s));
    mean4_DFA(s) = mean(DFA_tangential(win4_signal(1):win4_signal(2), s));
end

% correlation
[rho1, pval1] = corr(mean1_DFA(subj_in)', DFA_alpha_prestim_tang(subj_in)', 'type', 'Spearman'); % significant
[rho2, pval2] = corr(mean2_DFA(subj_in)', DFA_alpha_prestim_tang(subj_in)', 'type', 'Spearman'); % n.s.
[rho3, pval3] = corr(mean3_DFA(subj_in)', DFA_alpha_prestim_tang(subj_in)', 'type', 'Spearman'); % n.s.
[rho4, pval4] = corr(mean4_DFA(subj_in)', DFA_alpha_prestim_tang(subj_in)', 'type', 'Spearman'); % n.s.

% Bonferroni correction
pval_bonf = pval1*4;


%% control for SNR of SEP
load([savepath_pre 'SNR_all_subjects_MANUSCRIPT.mat']) % obtained from "N20_probe_manuscript_code_SNR_estimation.m"
load([savepath_pre 'DFA_timecourse_prestim_alpha_cca_MANUSCRIPT.mat'])

[rho, pval] = corr(SNR_mean_tang(subj_in), DFA_alpha_prestim_tang(subj_in)', 'type', 'Spearman'); % no correlation between SNR_SEP and DFA_alpha


    
%% relationship DFA exponents continuous alpha and DFA exponents early SEP
load('/data/pt_01972/Preproc_data/N20_study1/DFA_timecourse_cont_alpha_cca_MANUSCRIPT.mat')
load('/data/pt_01972/Preproc_data/N20_study1/DFA_exponents_with_permutations_and_CCA_averages_MANUSCRIPT.mat');

subj_in = [1:12, 14:17, 19:33]; % select subjects

% extract rms/ auc of DFA exponent time course
srate = 5000;
epoch_length_t = [-0.1 0.6];
win1_signal_ms = [20 25] - epoch_length_t(1)*1000; % in ms
win1_signal = win1_signal_ms * srate/1000;

win2_signal_ms = [25 30] - epoch_length_t(1)*1000; % in ms
win2_signal = win2_signal_ms * srate/1000;

win3_signal_ms = [30 35] - epoch_length_t(1)*1000; % in ms
win3_signal = win3_signal_ms * srate/1000;

win4_signal_ms = [35 40] - epoch_length_t(1)*1000; % in ms
win4_signal = win4_signal_ms * srate/1000;

mean1_DFA = []; mean2_DFA = []; mean3_DFA = []; mean4_DFA = [];
for s = 1:33    
    % average DFA exponents in time window of interest
    mean1_DFA(s) = mean(DFA_tangential(win1_signal(1):win1_signal(2), s));
    mean2_DFA(s) = mean(DFA_tangential(win2_signal(1):win2_signal(2), s));   
    mean3_DFA(s) = mean(DFA_tangential(win3_signal(1):win3_signal(2), s));
    mean4_DFA(s) = mean(DFA_tangential(win4_signal(1):win4_signal(2), s));
end

% correlation
[rho1, pval1] = corr(mean1_DFA(subj_in)', DFA_alpha_cont_tang(subj_in)', 'type', 'Spearman'); % significant
[rho2, pval2] = corr(mean2_DFA(subj_in)', DFA_alpha_cont_tang(subj_in)', 'type', 'Spearman'); % n.s.
[rho3, pval3] = corr(mean3_DFA(subj_in)', DFA_alpha_cont_tang(subj_in)', 'type', 'Spearman'); % n.s.
[rho4, pval4] = corr(mean4_DFA(subj_in)', DFA_alpha_cont_tang(subj_in)', 'type', 'Spearman'); % n.s.

    % Bonferroni correction
    pval_bonf = pval1*4;











