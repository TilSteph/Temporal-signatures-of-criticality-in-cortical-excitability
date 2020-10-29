%%% Temporal signatures of critical dynamics in cortical excitability as probed by early somatosensory responses %%%
% Authors: Tilman Stephani, Gunnar Waterstraat, Stefan Haufe, Garbiel Curio, Arno Villringer & Vadim V. Nikulin (2020)

%%% Code for model comparisons using ML-DFA (Ton & Daffertshofer, 2016) %%%


%% prepare environment
clear

addpath('/data/p_01972/EEGLAB/eeglab14_1_1b/');
addpath('/data/p_01972/Scripts/Code_data_manuscript_N20_power_law_dynamics/Functions/');
addpath(genpath('/data/p_01972/Other_toolboxes/FluctuationAnalysis_Daffertshofer/')) % Code published along Ton & Daffertshofer (2016)

eeglab


%% define names and paths
savepath_pre = '/data/pt_01972/Preproc_data/N20_study1/'; % first part of save directory
savepath_MANUSCRIPT = '/data/p_01972/Scripts/Code_data_manuscript_N20_power_law_dynamics/Data/'; % save path for preprocessed data that can be published according to the European GDPR

subj_dir = dir('/data/p_01972/RAW_EEG/N20_study1/*_*'); subj_names = {subj_dir.name};
subj_in = [1:12, 14:17, 19:33]; % select subjects
subj_names = subj_names(subj_in);
load([savepath_pre 'Classify_CCA_notask_30to200Hz_nonotch_stdzd_CP4.mat'])

srate = 5000;
epoch_length_ms = [-100, 600];
%intv_ms = [10, 50];
intv_ms = [-100, 200];
intv_pt = (intv_ms-epoch_length_ms(1))*srate/1000;
samples_sel = intv_pt(1)+1:intv_pt(2);
time_vec = intv_ms(1)+0.2:1/srate*1000:intv_ms(2);

Hurst_all = zeros(length(samples_sel),length(subj_names)); 
AIC_all = zeros(11,length(samples_sel),length(subj_names)); % model by sample by subject
BIC_all = zeros(11,length(samples_sel),length(subj_names));
for s = 1:length(subj_names)
    disp(subj_names{s})

    % load CCA data  
    subj = subj_names{s};
    pl_suffx = ['_notask_30to200Hz_nonotch_stdzd_CP4']; % according to loaded data
    savepath_full = [savepath_pre subj '/']; % complete save directory
    load([savepath_full 'CCA_' subj_names{s} pl_suffx '.mat'])

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

    % DFA parameters
    start = 7; 
    stop = 70;
    num_segment=30;
    start_fit = start; 
    stop_fit  = stop;
    fg=0; 
    %prune = [1, length(sig_x)-1]; % for comparison with DFA+

    fprintf(['\n' repmat('.',1,length(samples_sel)) '\n\n']);        
    parfor t = 1:length(samples_sel)
        fprintf('\b|\n'); % progress bar

        % extract signal for DFA
        sig_x = squeeze(CCA_comps(comps_tangential(s), samples_sel(t), :));

        % Daffertshofer DFA+
        [Hurst, details] = fluctuationAnalysis(cumsum(sig_x),[start, stop, num_segment],[],'DFA+');
        Hurst_all(t,s) = Hurst;
        AIC_all(:,t,s) = details.AICc;
        BIC_all(:,t,s) = details.BIC;
    end      
    fprintf('\n')
end

sig_x = squeeze(CCA_comps(comps_tangential(s), samples_sel(1), :));
[~, details] = fluctuationAnalysis(cumsum(sig_x),[start, stop, num_segment],[],'DFA+'); % to get example details 

% save
save([savepath_MANUSCRIPT 'DFA_Daffertshofer_TangCCA_10to50ms_MANUSCRIPT.mat'], 'Hurst_all', 'AIC_all', 'BIC_all', 'details', 'time_vec')



%% Evaluate DFA+ results
load([savepath_MANUSCRIPT 'DFA_Daffertshofer_TangCCA_10to50ms_MANUSCRIPT.mat'], 'Hurst_all', 'AIC_all', 'BIC_all', 'details', 'time_vec')

% plot model comparison over time course    
AIC_median = median(AIC_all,3);
BIC_median = median(BIC_all,3);
mod_labels = {'1+x' '1+x^2' '1+x+x^2' '1+x^3' '1+x+x^3' '1+x^2+x^3' '1+x+x^2+x^3' 'c(1)+c(2)*exp(c(3)*x)' ...
    'real((1/log(10))*log(c(1)*(1-exp(-c(2)*exp(log(10)*x)))+c(3)))' ...
    'piece-wise linear 2 sections' 'piece-wise linear 3 sections'};
cmap_lines = linspecer(length(mod_labels)); % colormap for model lines

% AICc
figure, plot(time_vec, AIC_median'), legend(mod_labels, 'location', 'southoutside') % -> model 1 and 2 best on average
    xlabel('time (ms)'), ylabel('median AICc'), title('ML-DFA model comparison (Daffertshofer) with empirical data (early SEPs)')
    set(gca, 'ColorOrder', cmap_lines);

% BIC
figure, plot(time_vec, BIC_median'), legend(mod_labels, 'location', 'southoutside') % -> model 1 and 2 best on average
    xlabel('time (ms)'), ylabel('median BIC'), title('ML-DFA model comparison (Daffertshofer) with empirical data (early SEPs)')
    set(gca, 'ColorOrder', cmap_lines);

% count aggregated model fits
%start_ix = find(time_vec==10); stop_ix = find(time_vec==50);
[~, mod_min] = min(BIC_median);
perc_loglog_median = sum(mod_min==1)/length(mod_min)*100; % 81.59%

% percentage retained data when excluding non log-log data
[~, mod_min] = min(BIC_all);
perc_loglog_all = sum(sum(mod_min==1))/numel(mod_min)*100; % 70.45%


% plot number linear DFA
figure, plot(time_vec, sum(~isnan(Hurst_all),2)) % -> minimum 17
xlabel('time (ms)'), ylabel('number of participants where log-log model fits best'), title('Frequency of ''log-log model outperforms other models'' in empirical data (early SEPs)')

% plot only exponents where DFA identified as linear
figure, plot(time_vec, nanmean(Hurst_all,2)) % -> comparable results (peaks at ~25 and 33 ms)
xlabel('time (ms)'), ylabel('Average DFA exponents (only included in average if log-log model accepted)'), title('DFA time course in early SEPs estimated with ML-DFA (Daffertshofer)')


    
    