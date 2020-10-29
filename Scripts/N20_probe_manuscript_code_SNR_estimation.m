%%% Temporal signatures of critical dynamics in cortical excitability as probed by early somatosensory responses %%%
% Authors: Tilman Stephani, Gunnar Waterstraat, Stefan Haufe, Garbiel Curio, Arno Villringer & Vadim V. Nikulin (2020)

%%% Code for SNR estimation of CCA components %%%


%% prepare environment
clear

addpath('/data/p_01972/EEGLAB/eeglab14_1_1b/');
addpath('/data/p_01972/Scripts/Code_data_manuscript_N20_power_law_dynamics/Functions/')

eeglab

% set eeglab options to double precision (important for CCA)
pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 1, 'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 0, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 0);


%% paths and names
savepath_pre = '/data/pt_01972/Preproc_data/N20_study1/'; % first part of save directory
datapath_pre = '/data/pt_01972/Preproc_data/N20_study1/'; % first part of data directory
savepath_figures = '/data/pt_01972/Results/figures/';

subj_dir = dir('/data/p_01972/RAW_EEG/N20_study1/*_*');
subj_names = {subj_dir.name};

% define conditions
cond = 'notask';

% load one EEG datasset and do CCA to get EEG.times etc
    EEG = pop_loadset('filename',['12_TS_notask_pchip_sr5kHz_30to200Hz_vi_averef_NEWnotch_ICA_removed.set'],'filepath', [datapath_pre '12_TS/']);
    EEG.data = double(EEG.data);    

    % CCA parameters
    epoch_length_t = [-0.1 0.6]; % in sec
    epoch_markers = {'A - Out', 'B - Out'};
    train_win = [5 80]; % "default": [5 80]
    eval_win = epoch_length_t*1000;
    num_comp = 4;
    flag_standardize = 1; % or 0; if 1, first two components' polarity is standardized
    
    % Do CCA
    [CCA_comps, Wst, Ast, R_cca, iCCA_tangential, iCCA_radial, EEG, ix_accepted] = apply_CCA(EEG, epoch_markers, train_win, eval_win, num_comp, flag_standardize); % do CCA

% windows for SNR estimation
win_signal_ms = [10 50] - epoch_length_t(1)*1000; % in ms
win_signal = win_signal_ms * EEG.srate/1000;

win_noise_ms = [-50 -10] - epoch_length_t(1)*1000; % in ms
win_noise = win_noise_ms * EEG.srate/1000;

% classification of components
load([savepath_pre 'Classify_CCA_notask_30to200Hz_nonotch_stdzd_CP4.mat'])



%% SNR in tangential CCA component and its relation to DFA exponents
SNR_mean_cca1 = zeros(length(subj_names),1);
SNR_mean_tang = zeros(length(subj_names),1);
for s = 1:length(subj_names)
    if strcmp(cond, 'task') && s==13 % subject 13 has no task data!
        continue
    end
    
    subj = subj_names{s};
    fprintf(['\n\n ------------ Subject ' subj '  condition: ' cond '   ------------ \n\n']);
   
    % load CCA data
    pl_suffx = ['_' cond '_30to200Hz_nonotch_stdzd']; % according to loaded data

    savepath_full = [savepath_pre subj '/']; % complete save directory
    load([savepath_full 'CCA_' subj_names{s} pl_suffx '.mat'])

    % estimate SNR
    rms_signal_cca1 = squeeze(rms(CCA_comps(1,win_signal(1):win_signal(2),:), 2)); 
    rms_signal_tang = squeeze(rms(CCA_comps(comps_tangential(s),win_signal(1):win_signal(2),:), 2)); % tangential CCA
        
    rms_noise_cca1 = squeeze(rms(CCA_comps(1,win_noise(1):win_noise(2),:), 2));
    rms_noise_tang = squeeze(rms(CCA_comps(comps_tangential(s),win_noise(1):win_noise(2),:), 2)); % tangential CCA
    
    SNR_tmp_cca1 = rms_signal_cca1./rms_noise_cca1;
    SNR_tmp_tang = rms_signal_tang./rms_noise_tang;
    
    SNR_mean_cca1(s) = mean(SNR_tmp_cca1);    
    SNR_mean_tang(s) = mean(SNR_tmp_tang);  
end

% save
save([datapath_pre 'SNR_all_subjects_MANUSCRIPT.mat'], 'SNR_mean_tang')
load([datapath_pre 'SNR_all_subjects_MANUSCRIPT.mat'], 'SNR_mean_tang')

% correlate with DFA exponents
load([savepath_pre 'DFA_exponents_with_permutations_and_CCA_averages_MANUSCRIPT.mat']); % DFA exponent time course of CCA components

subj_in = [1:12, 14:17, 19:33]; % select subjects
mean(SNR_mean_tang(subj_in))
std(SNR_mean_tang(subj_in))

% extraction windows
win_signal_ms = [10 50] - epoch_length_t(1)*1000; % in ms
win_signal = win_signal_ms * EEG.srate/1000;

% Average DFA exponents in time window
mean_DFA_tang = mean(DFA_tangential(win_signal(1):win_signal(2),subj_in), 1);

% Corr and plots
[rho, pval] = corr(mean_DFA_tang', SNR_mean_tang(subj_in), 'type', 'Spearman'); % significant;




%% SNR in thalamus activity
% thalamus-related components (30to200Hz_nonotch_CP4)
subj_P15 = [ 2, 4, 5, 8, 9, 11, 12, 21, 22, 24, 25, 26, 29]; % subject
comp_P15 = [ 4, 3, 3, 3, 4,  3,  3,  3,  4,  4,  4,  3,  3]; % CCA component with thalamus activity
pol_P15  = [-1, 1, 1, 1, 1,  1,  1, -1, -1,  1,  1,  1,  1]; % polarity of P15 peak

% windows for SNR estimation
win_signal_ms = [12 18] - epoch_length_t(1)*1000; % in ms
win_signal = win_signal_ms * EEG.srate/1000;

win_noise_ms = [-18 -12] - epoch_length_t(1)*1000; % in ms
win_noise = win_noise_ms * EEG.srate/1000;

% subject loop
SNR_mean_thal = zeros(length(subj_P15),1);
for s = 1:length(subj_P15)
    if strcmp(cond, 'task') && s==13 % subject 13 has no task data!
        continue
    end
    
    subj = subj_names{s};
    fprintf(['\n\n ------------ Subject ' subj '  condition: ' cond '   ------------ \n\n']);
   
    % load CCA data
    pl_suffx = ['_' cond '_30to200Hz_nonotch_stdzd_CP4']; % for this version thalamus was investigated
        
    savepath_full = [savepath_pre subj '/']; % complete save directory
    load([savepath_full 'CCA_' subj_names{s} pl_suffx '.mat'])

    % estimate SNR 
    rms_signal_thal = squeeze(rms(CCA_comps(comp_P15(s),win_signal(1):win_signal(2),:), 2)); % thalamic CCA
    rms_noise_thal = squeeze(rms(CCA_comps(comp_P15(s),win_noise(1):win_noise(2),:), 2)); % thalamic CCA
        
    SNR_tmp_thal = rms_signal_thal./rms_noise_thal;     
    SNR_mean_thal(s) = mean(SNR_tmp_thal);  
end

% descriptive statistics
mean(SNR_mean_thal)
std(SNR_mean_thal)
    


%% SNR of CNAP (median nerve)
% windows for SNR estimation
win_signal_ms = [5 8] - epoch_length_t(1)*1000; % in ms
win_signal = win_signal_ms * EEG.srate/1000;

win_noise_ms = [-8 -5] - epoch_length_t(1)*1000; % in ms
win_noise = win_noise_ms * EEG.srate/1000;

SNR_mean_CNAP = [];
for s = 1:length(subj_names)
    subj = subj_names{s};
    fprintf(['\n\n ------------ generate and save DFA time courses of subject ' subj '  condition: periphery ' cond '   ------------ \n\n']);
    
    savepath_full = [savepath_pre subj '/']; % complete save directory
    
    EEG = pop_loadset('filename',[subj_names{s} '_peri_notask_pchip_sr5kHz_HP70Hz_vi_notch.set'],'filepath', savepath_full);
        
    % Cut data into epochs
    epoch_length_t = [-0.1 0.6]; % in sec; shoudl be same length as where prune parameters are from
    [EEG, ix_accepted_events] = pop_epoch(EEG, {'A - Out', 'B - Out'}, epoch_length_t, 'newname', ['N20_probe_epochs'], 'epochinfo', 'yes'); % use second output later for choosing epochs (e.g. in behavioral data) AND prune parameters in DFA!!!!
        
    ichan_CNAP = 1;
    if s == 29, ichan_CNAP = 2; end
    
    % estimate SNR 
    rms_signal_CNAP = squeeze(rms(EEG.data(ichan_CNAP,win_signal(1):win_signal(2),:), 2)); % thalamic CCA
    rms_noise_CNAP = squeeze(rms(EEG.data(ichan_CNAP,win_noise(1):win_noise(2),:), 2)); % thalamic CCA
        
    SNR_tmp_CNAP = rms_signal_CNAP./rms_noise_CNAP;     
    SNR_mean_CNAP(s) = mean(SNR_tmp_CNAP);  
end

% descriptive statistics   
subj_in = [1:12, 14:17, 19:33];
mean(SNR_mean_CNAP(subj_in))
std(SNR_mean_CNAP(subj_in))
    




    
