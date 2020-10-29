%%% Temporal signatures of critical dynamics in cortical excitability as probed by early somatosensory responses %%%
% Authors: Tilman Stephani, Gunnar Waterstraat, Stefan Haufe, Garbiel Curio, Arno Villringer & Vadim V. Nikulin (2020)

%%% Code for figures in manuscript %%%
% Note: for export of the figures from Matlab, 'export_fig' was used (https://de.mathworks.com/matlabcentral/fileexchange/23629-export_fig)

%% prepare environment
addpath('/data/p_01972/EEGLAB/eeglab14_1_1b/');
addpath('/data/p_01972/Scripts/Misc/export_fig-master'); % contains export_fig library
addpath('/data/p_01972/Scripts/Code_data_manuscript_N20_power_law_dynamics/Functions/')

eeglab

% avoid "_" being interpreted as subscript
set(0, 'DefaulttextInterpreter', 'none')


%% paths and names
savepath_pre = '/data/pt_01972/Preproc_data/N20_study1/'; % first part of directory where processed data was saved
savepath_figures = '/data/pt_01972/Results/figures_manuscript/';
savepath_MANUSCRIPT = '/data/p_01972/Scripts/Code_data_manuscript_N20_power_law_dynamics/Data/'; % save path for preprocessed data that can be published according to the European GDPR

%% define colormaps
colormap_topo = diverging_map(linspace(0,1,100),[0.23 0.299 0.754], [0.706, 0.016, 0.15]);
colormap_sources = load('/data/p_01972/Scripts/Code_data_manuscript_N20_power_law_dynamics/Functions/cm17.mat', 'cm17a'); % colormap from grey to red


%% Fig 1: Ising model & paradigm
% The Ising model was simulated using the Metropolis algorithm; see function "metropolis_snap.m"

% Panel A: Ising model at different temperatures (=different system states)
load([savepath_MANUSCRIPT 'Ising_spin_snaps_256grid_1024snaps_Ts_low_crit_high_same_order.mat'])

h1 = figure;
for k = 1:3
    subplot(1,3,k), colormap(gray)
    imagesc(Ising(k).spin_snaps(:,:,end/2))
    set(gca,'visible','off'), set(gcf,'color','w')
    axis image
end

export_fig([figure_path 'Ising_3Temps_.0001_2.269_100_1024snaps.eps'], '-eps', '-transparent', h1)


% Panel B: five snapshots from the Ising model at the critical state
load([savepath_MANUSCRIPT 'Ising_spin_snaps_256grid_3072snaps_T2.2692.mat']) % contains "spin_snaps5"

n_snaps = size(spin_snaps5,3);
init_snaps = 1000; % 1000 instances to initialize dynamics
steps = floor((n_snaps-init_snaps)/4 * [0 1 2 3 4] +init_snaps ); 

spin_snaps_tmp = spin_snaps5;
h2 = figure;
for k = 1:5
    subplot(1,5,k), colormap(gray)
    imagesc(spin_snaps_tmp(:,:,steps(k)))
    set(gca,'visible','off'), set(gcf,'color','w')
    disp(steps(k))
    axis image
end

export_fig([figure_path 'Ising_5snaps_kTc_stepsize500.eps'], '-eps', '-transparent', h2)




%% Fig 2: Grand averages channel space and tangential CCA component + topographies + source reconstruction

% Panel A
% load data (SEP sensor space + tangential CCA)
load([savepath_MANUSCRIPT 'Sensor_grand_average_30to200Hz_nonotch.mat']) % Grand average in sensor space; also contains example single-subject EEG structure
load([savepath_MANUSCRIPT 'DFA_exponents_with_permutations_and_CCA_averages_MANUSCRIPT.mat']); % DFA exponent time courses (+permutations) and CCA averages
Ast_tangential(:,21) = -1 * Ast_tangential(:,21); Ast_radial(:,31) = -1 * Ast_radial(:,31); % adjust some polarities where automatic standardization did not work properly

% exclude 13 and 18
subj_in = [1:12, 14:17, 19:33]; % select subjects
mean_SEP_ex = mean_SEP_all(:,:,subj_in);

h2 = figure; hold on
plot(EEG.times, mean(mean_SEP_all(ismember({EEG.chanlocs.labels}, {'F4'}),:,subj_in), 3), 'b')
plot(EEG.times, mean(mean_SEP_all(ismember({EEG.chanlocs.labels}, {'CP4'}),:,subj_in), 3), 'g')
plot(EEG.times, mean(mean_SEP_all(ismember({EEG.chanlocs.labels}, {'P4'}),:,subj_in), 3), 'r')
ylim([-1.3 0.8])

yyaxis right
plot(EEG.times, mean(CCA_tangential_average(:,subj_in),2), 'k')
ylim([-1.3 0.8])
xlim([-50 150])
legend({'F4', 'CP4', 'P4', 'tangential CCA'}, 'location', 'southeast')

export_fig([savepath_figures 'Fig2_grand_average_sensor_F4_CP4_P4_tangCCA.eps'], '-eps', '-transparent', h2)

% Panel B
% plot sensor space topography
maplimits = [-0.7 0.7];
h4 = figure;
subplot(2,1,1)
topo_ms = 20;
topo_pt = (topo_ms - EEG.xmin*1000)*EEG.srate/1000;
topoplot(mean(mean_SEP_ex(:,topo_pt,:),3), EEG.chanlocs, 'colormap', colormap_topo, 'maplimits', maplimits);
colorbar
title([num2str(topo_ms) ' ms'])

% Panel C
% plot tangential CCA topography
subplot(2,1,2)
topoplot(mean(Ast_tangential(:,subj_in),2), EEG.chanlocs, 'colormap', colormap_topo, 'maplimits', maplimits);
colorbar,
title('tangential CCA')

export_fig([savepath_figures 'Fig2_topos_sensor_F4_CP4_P4_tangCCA.eps'], '-eps', '-transparent', '-painters', h4)


% Panels D & E
% plot result of source reconstruction: with Brainstorm GUI
% change some colorbar settings:
bst_colormaps('NewCustomColormap','source', [], colormap_sources) % change colormap
bst_colormaps('ConfigureColorbar', gcf, 'source', 'sloreta', []) % control colorbar ticks

print(gcf, [savepath_figures 'Sources_CCA_tangential_eLoreta_notask_without_13_and_18_0smooth_0amplitude_newcolormap.png'], '-dpng', '-r900'); % save plot
print(gcf, [savepath_figures 'Sources_CCA_tangential_eLoreta_notask_without_13_and_18_50smooth_95amplitude_newcolormap.png'], '-dpng', '-r900'); % save plot




%% Fig 3: Exemplary subject: single trials, SEP-over-trials time courses, log-log DFA plot, DFA time course

% Panel A
% Single trials of tangential CCA component of an exemplary subject
load([savepath_pre '27_TM/CCA_27_TM_notask_30to200Hz_nonotch_stdzd.mat']) % single-subject CCA data
x_lim = [-20 100];
k=2; % tangential CCA component
h = figure;
imagesc(time_vec, linspace(size(CCA_comps, 3), 1, size(CCA_comps, 3)), squeeze(CCA_comps(k,:,:))');
        colorbar; xlabel('time in ms'); ylabel('trials'); axis('xy')
        xlim(x_lim)
        caxis([-2 2])
        
export_fig([figure_path '27_TM_CCA_tangential_singletrials.eps'], '-eps', '-transparent', h)

% Panel B
% Plot SEP amplitude time courses
load([savepath_pre '27_TM/DFA_27_TM_notask_30to200Hz_nonotch_stdzd.mat']) % load single subject data containing DFA and amplitude information

h2 = figure;
subplot(3,1,1)
plot(amps_across_trials(1,:), 'b')
title(['at ' num2str(lats_ms(1)) ' ms'])

subplot(3,1,2)
plot(amps_across_trials(2,:), 'r')
title(['at ' num2str(lats_ms(2)) ' ms'])

subplot(3,1,3)
plot(amps_across_trials(3,:), 'k')
title(['at ' num2str(lats_ms(3)) ' ms'])

export_fig([savepath_figures 'Fig3b_amps_across_trials_20ms_25ms_29ms.eps'], '-eps', '-transparent', h2)


% Panel C
% plot log-log relationship of DFA
h3 = figure; hold on
plot(log10(time), log10(Alpha_lats(:,1)), '.b'); p1 = lsline;
plot(log10(time), log10(Alpha_lats(:,2)), '.r'); p2 = lsline;
plot(log10(time), log10(Alpha_lats(:,3)), '.k'); p3 = lsline;

xticks(log10([10 20 30 40 50 60 70]))
xticklabels({'10' '20' '30' '40' '50' '60' '70'})

yticks(log10(0.5:0.5:4))
yticklabels({0.5:0.5:4})

xlabel('window size tau (in trials)')
ylabel('fluctuation(tau)')
legend([p1 p2 p3], {'at 20 ms', 'at 25 ms', 'at 29 ms'}, 'location', 'southeast')

export_fig([savepath_figures 'Fig3c_loglog_20ms_25ms_29ms.eps'], '-eps', '-transparent', '-painters', h3)

% Panel D
% plot DFA exponents of all latencies (single subject)
ix_lats = find(ismember(round(EEG.times*1000)/1000, lats_ms));
sz = 10;
h4 = figure; hold on
plot(EEG.times, DFA_all_lats)
xlim([-20 100])

export_fig([savepath_figures 'Fig3d_DFA_time_course_27_TM.eps'], '-eps', '-transparent', h4)


% Panel E
% grand average DFA time course with SEP and cluster statistics
% load data
load([savepath_MANUSCRIPT 'DFA_exponents_with_permutations_and_CCA_averages_MANUSCRIPT.mat']); % DFA exponents
load([savepath_MANUSCRIPT 'Cluster_stats_tangentialCCA_without13and18_iter_1000_MANUSCRIPT.mat']); % results from cluster statistics

subj_in = [1:12, 14:17, 19:33]; % select subjects
epoch_length_ms = [-100 600];
time_vec = linspace(epoch_length_ms(1), epoch_length_ms(2), size(CCA_tangential_average,1));

x_lim = [-100 200];
y1_lim = [0.45 0.58];
y2_lim = [-1.5 1.5];

h5 = figure; hold on

% grey boxes to mark clusters
for idx = ix_survive_clusters
    box_start = time_vec(find(map_cluster==idx, 1, 'first'));
    box_end = time_vec(find(map_cluster==idx, 1, 'last'));
    p3 = fill([box_start,box_end,box_end,box_start],[y1_lim(1),y1_lim(1),y1_lim(2),y1_lim(2)], [.3 .3 .3], 'EdgeColor', 'none');
    %p3.FaceAlpha=0.3;
end

p1 = plot(time_vec, mean(DFA_tangential(:,subj_in),2));
ylabel('DFA exponent')
ylim(y1_lim)

yyaxis right
p2 = plot(time_vec, mean(CCA_tangential_average(:,subj_in),2));
xlim(x_lim), ylim(y2_lim)
xlabel('time (ms)'), ylabel('a.u.')
xticks([-100 -50 0 20 40 60 100 150 200])

legend([p1 p2 p3], {'DFA exponents', 'Average SEP (tangential CCA)', ['Significant clusters (ps<' num2str(p_mapthresh) ')']}, 'location', 'southoutside'), legend boxoff

export_fig([savepath_figures 'DFA_cluster_stats_longer.eps'], '-eps', '-transparent', '-painters', h5)



%% Fig 4: Control measures

% Panel A
% Single trials of CNAP (median nerve) of an exemplary subject
EEG_peri = pop_loadset('filename','12_TS_peri_notask_pchip_sr5kHz_HP70Hz_vi_notch.set','filepath', [savepath_pre '12_TS/']); % data from an exemplary subject
epoch_length_t = [-0.1 0.6]; % in sec
[EEG_peri, ix_accepted_events] = pop_epoch(EEG_peri, {'A - Out', 'B - Out'}, epoch_length_t, 'newname', ['N20_probe_epochs'], 'epochinfo', 'yes');
    
x_lim = [-20 100];
h = figure;
k=1; % CNAP channel
imagesc(EEG_peri.times, linspace(EEG_peri.trials, 1, EEG_peri.trials), squeeze(EEG_peri.data(k,:,:))');
        colorbar; xlabel('time in ms'); ylabel('trials'); axis('xy')
        xlim(x_lim)
        caxis([-10 10])
        
export_fig([savepath_figures '12_TS_CNAP_biceps_singletrials.eps'], '-eps', '-transparent', h)


% Panel B
% Grand average DFA exponents of CNAP + average SSEP
load([savepath_MANUSCRIPT 'DFA_exponents_time_CNAP_MANUSCRIPT.mat'])
subj_in = [1:12, 14:17, 19:33]; % select subjects

x_lim = [-100 200];
y_lim = [0.48 0.58];
h2 = figure;
p1 = plot(EEG.times, mean(DFA_exponents_time_peri(:, subj_in),2));
xlim(x_lim), xlabel('time (ms)')
ylim(y_lim), ylabel('DFA exponent')

yyaxis right
p2 = plot(EEG.times, mean(mean_SEP(:,subj_in),2));
ylabel('Amplitude (µV)')
ylim([-1.3 1.7]), 

legend([p1 p2], {'DFA exponents', 'CNAP'}, 'location', 'northeast'), legend boxoff
set(gcf,'color','white') % white background

export_fig([savepath_figures 'DFA_timecourse_biceps_plus_CNAP.eps'], '-eps', '-transparent', h2)


% Panel C
% Single trials thalamus CCA component (examplary subject)
load([savepath_pre '12_TS/CCA_12_TS_notask_30to200Hz_nonotch_stdzd.mat']) % exemplary subject
thal_comp_12 = 3;
x_lim = [-20 100];
h3 = figure;
imagesc(time_vec, linspace(size(CCA_comps, 3), 1, size(CCA_comps, 3)), -1*squeeze(CCA_comps(thal_comp_12,:,:))'); % reverse arbitrary polarity for illustration purposes
        colorbar; xlabel('time in ms'); ylabel('trials'); axis('xy')
        xlim(x_lim)
        caxis([-2 2])
        
export_fig([savepath_figures '12_TS_CCA_thalamus_singletrials_revpol.eps'], '-eps', '-transparent', h3)

        
% Panel D
% Grand average DFA exponents of thalamus-related CCA components + SEP
load([savepath_MANUSCRIPT 'DFA_exponents_with_permutations_and_CCA_averages_MANUSCRIPT.mat']);
load([savepath_MANUSCRIPT 'DFA_exponents_time_thalamus_MANUSCRIPT.mat'])

x_lim = [-100 200];
y_lim = [0.48 0.58];

h4 = figure;
p1 = plot(EEG.times, mean(DFA_cca_thal,2));
xlim(x_lim), xlabel('time (ms)')
ylim(y_lim), ylabel('DFA exponent')

yyaxis right
p2 = plot(EEG.times, mean(mean_CCA_thalamus,2));
ylabel('Amplitude (a.u.)')

legend([p1 p2], {'DFA exponents', 'P15 component'}, 'location', 'northeast'), legend boxoff
set(gcf,'color','white') % white background

export_fig([savepath_figures 'DFA_timecourse_thalamus_plus_SEP.eps'], '-eps', '-transparent', h4)


% Panel E
% Average spatial pattern of thalamus-related CCA components
Ast_thal = [];
for s = 1:length(subj_P15)
    Ast_thal(:,s) = Ast_all(comp_P15(s), :, subj_P15(s)) * pol_P15(s); % and standardize polarity
end

h5=figure;
maplimits = [-0.7 0.7];
topoplot(mean(Ast_thal,2), EEG.chanlocs, 'colormap', colormap_topo, 'maplimits', maplimits);
colorbar

export_fig([savepath_figures 'Spatial_pattern_thalamus_CCA_average.eps'], '-eps', '-painters', '-transparent', h5)


% Panel F
% SNR simulation
load([savepath_MANUSCRIPT 'DFA_mixing_distribution.mat'])

% only display the relevant parameter regions of the simulations
zoom_SNR = 500;
zoom_DFA = 30;
h6 = figure;
imagesc(mean(SNR_test(1:zoom_SNR,:),2), mean(exponent_true_emp(1:zoom_DFA,:),2), mean(exponent_mix(1:zoom_SNR, 1:zoom_DFA,:),3)')
axis xy

% mark coordinates of empirical SNR and DFA exponent
DFA_emp = 0.575; SNR_emp = 1.64;
SNR_all = mean(SNR_test(1:zoom_SNR,:),2);
DFA_mixed_zoomed = mean(exponent_mix(1:zoom_SNR, 1:zoom_DFA,:),3);
DFA_true_zoomed = mean(exponent_true_emp(1:zoom_DFA,:),2);
[~, ix_SNR] = min(abs(SNR_all-SNR_emp));
[~, ix_mixed] = min(abs(DFA_mixed_zoomed(ix_SNR,:)-DFA_emp));

hold on
plot(SNR_emp,  DFA_true_zoomed(ix_mixed), 'or')

% add labels
xlabel('SNR')
ylabel('true DFA exponent')
h = colorbar; set(get(h,'label'),'string','DFA exponent of mixed signal');
set(gcf,'color','white') % white background

export_fig([savepath_figures 'SNR_mixing_simulation.eps'], '-eps', '-transparent', h6)




%% Fig. 5: Relationship between pre-stimulus alpha and early SEP

% Panel A: Amplitude relation between pre-stimulus alpha and SEP (both tangential CCA)
load([savepath_MANUSCRIPT 'Predict_N20_from_prestim_alpha_200to10ms_CCA_filter_MANUSCRIPT.mat'])

subj_in = [1:12, 14:17, 19:33]; % select subjects
TangCCA_bins_mean = TangCCA_bins_mean(:,:,subj_in); 

epoch_length_ms = [-100 600];
time_vec = linspace(epoch_length_ms(1), epoch_length_ms(2), size(TangCCA_bins_mean,1));

h5a = figure;
subplot(1,2,1), hold on
p1 = plot(time_vec, mean(TangCCA_bins_mean(:,1,:),3));
p2 = plot(time_vec, mean(TangCCA_bins_mean(:,2,:),3));
p3 = plot(time_vec, mean(TangCCA_bins_mean(:,3,:),3));
p4 = plot(time_vec, mean(TangCCA_bins_mean(:,4,:),3));
p5 = plot(time_vec, mean(TangCCA_bins_mean(:,5,:),3));
legend([p1 p2 p3 p4 p5], {'bin1', 'bin2', 'bin3', 'bin4', 'bin5'}, 'location', 'southoutside')
xlim([-20 100])

% zoom in
subplot(1,2,2), hold on
p1 = plot(time_vec, mean(TangCCA_bins_mean(:,1,:),3));
p2 = plot(time_vec, mean(TangCCA_bins_mean(:,2,:),3));
p3 = plot(time_vec, mean(TangCCA_bins_mean(:,3,:),3));
p4 = plot(time_vec, mean(TangCCA_bins_mean(:,4,:),3));
p5 = plot(time_vec, mean(TangCCA_bins_mean(:,5,:),3));
legend([p1 p2 p3 p4 p5], {'bin1', 'bin2', 'bin3', 'bin4', 'bin5'}, 'location', 'southoutside')
xlim([18.5 22])

export_fig([savepath_figures 'Fig9a_relation_DFA_prestim_alpha_and_SEP_bins.eps'], '-eps', '-transparent', '-painters', '-nocrop', h5a)


% Panel B: DFA relation between pre-stimulus alpha and SEP (both tangential CCA)
load([savepath_MANUSCRIPT 'DFA_timecourse_prestim_alpha_cca_MANUSCRIPT.mat'])
load([savepath_MANUSCRIPT 'DFA_exponents_with_permutations_and_CCA_averages_MANUSCRIPT.mat']);

subj_in = [1:12, 14:17, 19:33]; % select subjects
srate = 5000; % sampling rate

% extract rms/ auc of DFA exponent time course
epoch_length_t = [-0.1 0.6];
win1_signal_ms = [20 25] - epoch_length_t(1)*1000; % in ms
win1_signal = win1_signal_ms * srate/1000;

mean1_DFA = []; %rms1_DFA = []; 
for s = 1:33    
    %rms1_DFA(s) = rms(DFA_tangential(win1_signal(1):win1_signal(2), s));
    mean1_DFA(s) = mean(DFA_tangential(win1_signal(1):win1_signal(2), s)); % changed to mean(); TS 04/2020
end

h5b = figure; hold on
x_lim = [0.45 0.8];
y_lim = [0.47 0.8];
%scatter(DFA_alpha_prestim_tang(subj_in), rms1_DFA(subj_in), '.')
scatter(DFA_alpha_prestim_tang(subj_in), mean1_DFA(subj_in), '.') % changed to mean(); TS 04/2020
xlim(x_lim)
ylim(y_lim)
%[B1,BINT1,R1] = regress(rms1_DFA(subj_in)', [ones(size(DFA_alpha_prestim_tang(subj_in)')), DFA_alpha_prestim_tang(subj_in)']);
[B1,BINT1,R1] = regress(mean1_DFA(subj_in)', [ones(size(DFA_alpha_prestim_tang(subj_in)')), DFA_alpha_prestim_tang(subj_in)']);
plot(DFA_alpha_prestim_tang(subj_in), B1(1)+B1(2)*DFA_alpha_prestim_tang(subj_in))
xlabel('DFA exponents of prestimulus alpha')
ylabel('average DFA exponents of SEP')

export_fig([savepath_figures 'Fig5b_relation_DFA_prestim_alpha_and_SEP_mean.eps'], '-eps', '-transparent', '-nocrop', h5b)






