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
savepath_MANUSCRIPT = '/data/p_01972/Scripts/Code_data_manuscript_N20_power_law_dynamics/Data/'; % save path for preprocessed data that can be published according to the European GDPR

subj_dir = dir('/data/p_01972/RAW_EEG/N20_study1/*_*');
subj_names = {subj_dir.name};

% define conditions
cond = 'notask';
    

%% subject loop
nreps = 1000; % repetitions for DFA surrogates
DFA_alpha_cont_tang_perm = zeros(nreps,length(subj_names));
DFA_alpha_cont_tang = []; DFA_alpha_cont_rad = [];
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
         
    % load previously calculated CCA data
    load([savepath_full 'CCA_' subj_names{s} '_' cond '_30to200Hz_nonotch_stdzd_CP4.mat'])    
    load([savepath_pre 'Classify_CCA_notask_30to200Hz_nonotch_stdzd_CP4.mat']) % classification of components
    
    
    %% trial loop
    filt_hz = [8 13]; % [4 7]; [8 10]; [10 13]; [15 20]; [20 30];
    [b,a] = butter(2, filt_hz/(EEG.srate/2));
    
    sig = EEG.data(:,:);
    sig = filtfilt(b, a, double(sig)'); % filter in alpha range     
    
    % apply CCA weights to prestimulus alpha
    sig_cca = Wst(:,[comps_tangential(s) comps_radial(s)])' * sig'; % only look at tangential and radial CCA components

    % get envelope of band-pass filtered data
    sig_cca_env = abs(hilbert(sig_cca'))'; 
       
    EEG.data = sig_cca_env; % put filtered signal (CCA components) back in EEG structure
    EEG = eeg_checkset(EEG); % correct other EEG fields    
    
    
    %% DFA
    start_sec = 5; % corresponding to ~7 trials
    stop_sec = 50; % corresponding to ~70 trials
    
    % DFA parameters
    start = start_sec*EEG.srate;
    stop = stop_sec*EEG.srate;
    num_segment = 30;
    start_fit = start;
    stop_fit  = stop;
    fg=0;   
    
    % prune paramters (account for breaks in data)
    boundary_lats = [EEG.event(find(ismember({EEG.event.type}, 'boundary'))).latency];
    segment_pnts = [1, ceil(boundary_lats), EEG.pnts]; % round up
    segment_dur = diff(segment_pnts)-1; % durations until one sample before boundary element
    prune = [segment_pnts(1:end-1); segment_dur]';  
    
    [DFA_alpha_cont_tang(s),Amplitude,Alpha,time,st,epochs]=dfa_2018(EEG.data(1,:),start,stop,num_segment,start_fit,stop_fit,prune, fg); % for tangential sources
    [DFA_alpha_cont_rad(s),Amplitude,Alpha,time,st,epochs]=dfa_2018(EEG.data(2,:),start,stop,num_segment,start_fit,stop_fit,prune, fg); % for radial sources    
    
    % permuted tangential
    parfor k = 1:nreps
        sig_perm = EEG.data(1,randperm(EEG.pnts)); % shuffle samples
        [DFA_alpha_cont_tang_perm(k,s), Amplitude,Alpha,time,st,epochs]=dfa_2018(sig_perm,start,stop,num_segment,start_fit,stop_fit,prune, fg); % with prune parameters
    end
    
end

% save
save([savepath_MANUSCRIPT 'DFA_timecourse_cont_alpha_cca_MANUSCRIPT.mat'], 'DFA_alpha_cont_tang', 'DFA_alpha_cont_rad', 'DFA_alpha_cont_tang_perm')


% relationship to prestimulus alpha DFA
load([savepath_MANUSCRIPT 'DFA_timecourse_prestim_alpha_cca_MANUSCRIPT.mat']) % contains 'DFA_alpha_prestim_tang'
[r, pval] = corr(DFA_alpha_cont_tang(subj_in)', DFA_alpha_prestim_tang(subj_in)', 'type', 'Spearman'); % r ~ 0.9



%% Surrogate test for DFA exponents of continuous alpha
n_reps = 1000; % repetitions for DFA surrogates
n_data = 3600000; % ~12 min
srate = 5000;

DFA_exponent_sim = zeros(n_reps,1);
fprintf('\n');
parfor k = 1:n_reps
    fprintf('\b.\n');
    sig_raw = rand(1, n_data);
    
    % filter signal and get envelope
    filt_hz = [8 13];
    [b,a] = butter(2, filt_hz/(EEG.srate/2));
    
    sig_filt = filtfilt(b, a, sig_raw); % filter in alpha range 
    sig_env = abs(hilbert(sig_filt)); % envelope
    
    % DFA
    start_sec = 5; % corresponding to ~7 trials
    stop_sec = 50; % corresponding to ~70 trials
    
    start = start_sec*srate;
    stop = stop_sec*srate;
    num_segment = 30;
    start_fit = start;
    stop_fit  = stop;
    prune = [1 n_data-1]; % assume no breaks
    fg=0;   
    
    [DFA_exponent_sim(k),Amplitude,Alpha,time,st,epochs]=dfa_2018(sig_env, start,stop,num_segment,start_fit,stop_fit,prune, fg);
    
end
fprintf('\n\n');

save([savepath_MANUSCRIPT 'DFA_surrogates_for_cont_alpha_MANUSCRIPT.mat'], 'DFA_exponent_sim', 'n_data', 'n_reps')

% look at distribution
figure, hist(DFA_exponent_sim)

% descriptive statistics
load([savepath_MANUSCRIPT 'DFA_timecourse_cont_alpha_cca_MANUSCRIPT.mat'], 'DFA_alpha_cont_tang', 'DFA_alpha_cont_rad', 'DFA_alpha_cont_tang_perm')
subj_in = [1:12,14:17,19:33]; 
mean(DFA_alpha_cont_tang(subj_in))
mean(DFA_exponent_sim)

% statistical significance
[H,P,CI,STATS] = ttest(DFA_alpha_cont_tang(subj_in)', mean(DFA_exponent_sim));  

