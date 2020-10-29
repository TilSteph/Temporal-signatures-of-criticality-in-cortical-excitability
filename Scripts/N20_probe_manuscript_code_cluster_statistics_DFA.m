%%% Temporal signatures of critical dynamics in cortical excitability as probed by early somatosensory responses %%%
% Authors: Tilman Stephani, Gunnar Waterstraat, Stefan Haufe, Garbiel Curio, Arno Villringer & Vadim V. Nikulin (2020)

%%% Code for cluster statistics of DFA exponent time courses %%%

% Procedure:
% Compare clusters in t value time course of "emp. DFA vs average surrogate
% data" to clusters in t value time courses of null distribution "one
% surrogate vs surrogate average" (we assume the average surrogate data to 
% be the "true" DFA exponents at stochastical independence over trials.


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


%% Load DFA exponent time courses
load([savepath_MANUSCRIPT 'DFA_exponents_with_permutations_and_CCA_averages_MANUSCRIPT.mat']) % contains DFA exponent time courses plus surrogate DFA exponents (by permuting the trial order 1000 times)
subj_in = [1:12, 14:17, 19:33]; % select subjects

DFA_exponents_time = DFA_tangential(:,subj_in);
DFA_exponents_time_perm = DFA_tangential_perm(:,subj_in,:);


%% Identify clusters in empirical data
nreps = size(DFA_exponents_time_perm,3);
p_prethresh = 0.001; % in order to derive clusters
p_mapthresh = 0.001; % alpha level of significance of clusters

% pre-threshold for clusters
H = zeros(1,size(DFA_exponents_time,1));
T = zeros(1,size(DFA_exponents_time,1));
parfor i = 1:size(DFA_exponents_time,1) % loop through all empirical samples
    [H(i), P(i), ~, Stats] = ttest(DFA_exponents_time(i,:)', squeeze(mean(DFA_exponents_time_perm(i,:,:),3))', 'alpha', p_prethresh); % compare with mean of DFA time series at random; two-sample t-test
    T(i) = Stats.tstat;
end

% find clusters of significant T-values (according to pre-threshold)
[map_cluster, n_cluster] = bwlabel(H);

% Sum of samples and T-values of each cluster
clustcount = zeros(1,n_cluster);
T_sum   = zeros(1,n_cluster);
for i=1:n_cluster
    clustcount(i) = sum(map_cluster(:)==i); % samples per cluster
    T_sum(i)   = sum(T(map_cluster(:)==i)); % T-value sum
end
       

%% Generate null distribution      
DFA_exponents_time_perm_MEAN = squeeze(mean(DFA_exponents_time_perm(:,:,:),3));

max_T_sum = zeros(1,nreps);

parfor k = 1:nreps            
    DFA_exponents_time_perm_tmp = DFA_exponents_time_perm(:,:,k);            

    % pre-threshold for clusters
    H_null = zeros(size(DFA_exponents_time,1), 1);
    T_null = zeros(size(DFA_exponents_time,1), 1);
    for i = 1:size(DFA_exponents_time,1) % loop through all empirical samples
        [H_null(i), ~, ~, Stats] = ttest(DFA_exponents_time_perm_tmp(i,:)', DFA_exponents_time_perm_MEAN(i,:)', 'alpha', p_prethresh); % compare with mean of DFA time series at random; two-sample t-test
        T_null(i) = Stats.tstat;
    end

    % find clusters of significant T-values (according to pre-threshold)
    [map_cluster_null, n_cluster] = bwlabel(H_null);

    % Sum of samples and T-values of each cluster
    clustcount = zeros(1,n_cluster);
    T_sum_null = zeros(1,n_cluster);
    if n_cluster>=1
        for n = 1:n_cluster
            clustcount(n) = sum(map_cluster_null(:) == n); % samples per cluster
            T_sum_null(n) = sum(T_null(map_cluster_null(:) == n)); % T-value sum
        end

        max_T_sum(k) = max(T_sum_null); % find maximum t value cluster
    
    end

    if mod(k, 50) == 0, fprintf('permutation #%d \n', k), end % communicate with the operator       
end

% look at distribution of largest cluster t values
figure, hist(max_T_sum, 50)

% set T-value of map threshold
max_T_sum_sort = sort(max_T_sum(1,:)); % use null distribution from largest cluster t values
T_crit = max_T_sum_sort(round((1-p_mapthresh)*length(max_T_sum_sort))); % one-sided test (only for high T-values)

% "Significant" clusters (DFA exponents higher than H_random)
ix_survive_clusters = find(T_sum > T_crit);
ix_survive_samples = find(ismember(map_cluster, ix_survive_clusters));

% Mark surviving samples in DFA time course
mean_DFA_time = mean(DFA_exponents_time,2);
time_vec = linspace(-100, 600, 3500);
figure; hold on
plot(time_vec, mean_DFA_time)
plot(time_vec(ix_survive_samples), mean_DFA_time(ix_survive_samples), '.r')
plot(time_vec, squeeze(mean(mean(DFA_exponents_time_perm(:,:,:),3),2)), 'k')
xlim([-20 100])

% save cluster stats
save([savepath_MANUSCRIPT 'Cluster_stats_tangentialCCA_without13and18_iter_1000_MANUSCRIPT.mat'], 'ix_survive_samples', 'ix_survive_clusters', 'T_sum', 'T_crit', 'map_cluster', 'p_prethresh', 'p_mapthresh');
 