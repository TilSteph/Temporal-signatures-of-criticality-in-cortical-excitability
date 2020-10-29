%%% Temporal signatures of critical dynamics in cortical excitability as probed by early somatosensory responses %%%
% Authors: Tilman Stephani, Gunnar Waterstraat, Stefan Haufe, Garbiel Curio, Arno Villringer & Vadim V. Nikulin (2020)

%%% Code for data preprocessing %%%

% Note that the dataset which is preprocessed here consisted of different
% paradigms of which only the first one is addressed in the present
% manuscript ("Passive stimulation sequence").
% 1) Passive stimulation sequence
% 2) Somatosensory discrimination task (reported elsewhere: Stephani et
% al., in prepration)


%% environment preparation
clear

addpath('/data/p_01972/EEGLAB/eeglab14_1_1b/');
addpath('/data/p_01972/Scripts/Code_data_manuscript_N20_power_law_dynamics/Functions/')

eeglab
% set options: single precision for preprocessing
pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 1, 'option_single', 1, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 0, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 0);


%% Define paths and names
savepath_pre = '/data/pt_01972/Preproc_data/N20_study1/'; % first part of data directory
datapath_raw = '/data/p_01972/RAW_EEG/N20_study1/'; % first part of save directory

subj_dir = dir('/data/p_01972/RAW_EEG/N20_study1/*_*');
subj_names = {subj_dir.name};

ses_names = {};
% get ses file names
for s=1:length(subj_names)
    data_path = [datapath_raw subj_names{s} '/'];
    ses_files_subj = dir([data_path '*.ses']);
    ses_files_subj_name = sort({ses_files_subj.name});
    if s<=8
        ses_names{s} = ses_files_subj_name{1};
    else % there is also resting-state data
        ses_names{s} = ses_files_subj_name{2};
    end
end

electr_pos_file = '/data/p_01972/EEGLAB/standard-10-5-cap385_plusEOG.elp'; % adapted electrode positions for plotting


%% subject loop
for s = 1:33
    subj = subj_names{s};
    fprintf(['\n\n ------------   SUBJECT ' subj '   ------------ \n\n']);
    
    savepath_full = [savepath_pre subj '/']; % complete save directory
    
    %% Load data
    % Define path
    data_path = [datapath_raw subj '/']; % ses file location
    
    subj = [subj '_pchip'];
    
    % .ses file
    ses_file = ses_names{s};
    
    % Imports
    if strcmp(subj_names{s}, '12_TS') % recording was splitted for 12_TS 
        ses_file_notask = 'NeurOne-2018-08-03T120052.ses';
        ses_file_task = 'NeurOne-2018-08-03T121851.ses';
        EEG1 = pop_readneurone([data_path ses_file_notask], 1, '');
        EEG2 = pop_readneurone([data_path ses_file_task], 1, '');
        EEG = pop_mergeset(EEG1, EEG2);
        clear EEG1 EEG2
    else
        EEG = pop_readneurone([data_path ses_file], 1, '');
    end    
    EEG.etc = []; % clear this info (contains NeurOne infos); not necessary and interferes with following eeglab functions!

    
    %% Find event latencies
    % where to cut (triggers) 
    ix_A_Out = find(ismember({EEG.event.type}, 'A - Out')); % stimulation events A
    lat_A_Out = [EEG.event(ix_A_Out).latency];

    ix_B_Out = find(ismember({EEG.event.type}, 'B - Out')); % Stimulation B only in 2nd condition
    lat_B_Out = [EEG.event(ix_B_Out).latency];

    lat_all_Out = sort([lat_A_Out, lat_B_Out]);
    
        
    %% Interpolation of stimulation artifact using pchip = "cubic monotonous hermite spline interpolation"      
    % cut out samples (from -2 to 4 ms)
    sr = EEG.srate;
    t_cut = [-2 4]; % interpolated interval around triggers; in ms 
    pt_cut = ceil(t_cut/1000*sr); % in sampling points

    % Piecewise Cubic Hermite Interpolating Polynomial (PCHIP) + replace EEG data
    n_pt_fit = 5; %+1;  number of samples before and after cut used for interpolation fit
    
    x_fit_raw = [pt_cut(1)-n_pt_fit : 1 : pt_cut(1), pt_cut(2) : 1 : pt_cut(2)+n_pt_fit];
    x_sr_raw = [pt_cut(1) : 1 : pt_cut(2)]; % points to be interpolated; in pt

    for i = 1:length(lat_all_Out) % loop through all stimulation events
        x_fit = lat_all_Out(i) + x_fit_raw; % fit point latencies for this event
        x_sr = lat_all_Out(i) + x_sr_raw; % latencies for to-be-interpolated data points

        for c = 1:size(EEG.data,1) % loop through all channels
            y_fit = EEG.data(c, x_fit); % y values to be fitted
            %y_interp = pchip(x_fit, y_fit, EEG.data(c, x_sr-x_sr_raw)); % calculate pchip, obtain values in t_cut interval
            y_interp = pchip(x_fit, y_fit, x_sr); % calculate pchip, obtain values in t_cut interval
            EEG.data(c, x_sr) = y_interp; % replace in EEG data
        end

        if mod(i, 100) == 0 % show message every 100th trial
            fprintf('stimulation event %d \n', i)
        end
    end


    %% Take only experimental portions of the data (no beginning, end, or breaks)
    % (to save time during visual inspection)
    
    start_notask_tr = 11;
    end_notask_tr = 12;
    start_task_tr = 31;
    end_task_tr = 32;
    
    start_notask_lat = [EEG.event(  find(ismember({EEG.event.type}, num2str(start_notask_tr)), 1, 'last')  ).latency]; % get latency
    end_notask_lat = [EEG.event(  find(ismember({EEG.event.type}, num2str(end_notask_tr)), 1, 'last')  ).latency]; % get latency    
    start_task_lat = [EEG.event(  find(ismember({EEG.event.type}, num2str(start_task_tr)), 1, 'last')  ).latency]; % get latency
    end_task_lat = [EEG.event(  find(ismember({EEG.event.type}, num2str(end_task_tr)), 1, 'last')  ).latency]; % get latency
        
    offset = 3; % in sec
    if ~strcmp(subj_names{s}, '13_FB')
        EEG = pop_select(EEG, 'point', [start_notask_lat - offset*EEG.srate, end_notask_lat + offset*EEG.srate; start_task_lat - offset*EEG.srate, end_task_lat + offset*EEG.srate]); % take only experimental data    
    elseif strcmp(subj_names{s}, '13_FB')
        EEG = pop_select(EEG, 'point', [start_notask_lat - offset*EEG.srate, end_notask_lat + offset*EEG.srate]); % only no task condition
    end
    
    % update event latencies
    ix_A_Out = find(ismember({EEG.event.type}, 'A - Out')); % stimulation events A
    lat_A_Out = [EEG.event(ix_A_Out).latency];

    ix_B_Out = find(ismember({EEG.event.type}, 'B - Out')); % Stimulation B only in 2nd condition
    lat_B_Out = [EEG.event(ix_B_Out).latency];
    
    lat_all_Out = sort([lat_A_Out, lat_B_Out]);
    
    
    % exclude breaks (assume notask is always before task condition)
    notask_break_1_start = lat_all_Out(500);
    notask_break_1_end = lat_all_Out(501);
    
    if ~strcmp(subj_names{s}, '13_FB')
        task_break_1_start = lat_all_Out(1200);
        task_break_1_end = lat_all_Out(1201);

        task_break_2_start = lat_all_Out(1400);
        task_break_2_end = lat_all_Out(1401);

        task_break_3_start = lat_all_Out(1600);
        task_break_3_end = lat_all_Out(1601);

        task_break_4_start = lat_all_Out(1800);
        task_break_4_end = lat_all_Out(1801);
    end
    
    if ~strcmp(subj_names{s}, '13_FB')
    break_segments = [notask_break_1_start, notask_break_1_end; ...
                    task_break_1_start, task_break_1_end; ...
                    task_break_2_start, task_break_2_end; ...
                    task_break_3_start, task_break_3_end; ...
                    task_break_4_start, task_break_4_end];    
    elseif strcmp(subj_names{s}, '13_FB') % only notask condition
        break_segments = [notask_break_1_start, notask_break_1_end];
    end
    
    offset = 3; % in sec
    for q=1:size(break_segments, 1); if break_segments(q,2)-offset*EEG.srate < break_segments(q,1)+offset*EEG.srate; fprintf(['break shorter than chosen offset (break ' num2str(q) ')!\n']); end; end % check offsets in breaks
    
    if ~strcmp(subj_names{s}, '01_BT') % 01_BT: no breaks
        EEG = pop_select(EEG, 'nopoint', break_segments + [offset*EEG.srate, -offset*EEG.srate]); % take only experimental data
    end
        
    
    %% Exclude peripheral channels
    indxPeri = find(ismember({EEG.chanlocs.labels}, {'Thumb' 'Biceps'})); % indices of EOG channels
    EEG_peri = pop_select(EEG,'channel',  indxPeri); % only peripheral channels
    EEG = pop_select(EEG,'nochannel',  indxPeri);
    
    
    %% Optional: bad channel interpolation (if we already know those)
    if s == 1 %subject 01
        chanlocs_before = {EEG.chanlocs.labels};
        bad_elec = find(ismember({EEG.chanlocs.labels}, 'FC3'));
        EEG.etc.channel_interpolation.interpolated_electrodes = {EEG.chanlocs(bad_elec).labels}; % remember which electrode was interpolated
        
        [EEG, EEG.etc.channel_interpolation.com] = pop_interp(EEG, bad_elec, 'spherical');
        chanlocs_after = {EEG.chanlocs.labels}; % matches chanlocs_before!
    end
    
    EEG_N20 = EEG; % maintain/ copy for later analyses

    %% Down-sampling (helps ICA)
    down_sr = 250;
    EEG = pop_resample(EEG, down_sr); % down-sample to 250 Hz (anti-aliasing high-pass filter included)

    %% Visual data inspection (filter applied only for detection; prune parameters applied to original data)
    % Filter for visual inspection
    [b,a] = butter(2,[1 45]/(EEG.srate/2)); % as in Elena?s script
    EEG.data = filtfilt(b, a, double(EEG.data)')'; % dataset for artifact detection and later ICA
    
    % notch filter (only for visual inspection)
    bsFilter = [49 51];
    [b_notch, a_notch] = butter(2, bsFilter/(EEG.srate/2),'stop');
    EEG.data = filtfilt(b_notch, a_notch, double(EEG.data)')'; % dataset for artifact detection and later ICA    

    % A) Exclude (really) bad channels
        figure; pop_spectopo(EEG, 1, [], 'EEG' , 'percent', 100, 'freqrange',[1 100],'electrodes','off'); % spectrogram

        % pop_eegplot( EEG, 1, 1, 1); % show time series

        indxEog = find(ismember({EEG.chanlocs.labels}, {'VEOG1', 'VEOG2', 'HEOG1', 'HEOG2'})); % indices of EOG channels
        fprintf('\n  Do not remove EOG Channel (index:  %d)   \n', indxEog)

        EEG = pop_select(EEG); %( EEG,'nochannel',{'CP1'});
        reply = input('Remove channels, and band segments, then type c to continue: ','s');

        % remove same channels in N20 data (for later analyses)
        EEG_N20 = pop_select(EEG_N20); % (EEG, 'nochannel',{'FDI'});
        reply = input('Remove same channels in dataset for actual N20 analyses, then type c to continue: ','s');

        fprintf('  Number of left channels: %d   \n',EEG.nbchan)


    % B) Mark bad segments but do not exclude them (only mark)
        % set some plot properties
        marked_segments = []; % initialize variable for marks
        eegplotoptions = { 'events', EEG.event, 'eloc_file', EEG.chanlocs};

        % specify command for 'REJECT' button: save start points and durations of selected segments
        command = ['marked_segments.seg_start = TMPREJ(:,1);' ...
                   'marked_segments.seg_dur = TMPREJ(:,2) - TMPREJ(:,1);']; % TMPREJ = [start end R G B channels...]

        % Plot data and change the label of the 'REJECT' button to 'mark segments'
        eegplot( EEG.data, 'srate', EEG.srate, 'title', 'Scroll channel activities -- eegplot()', ...
                      'limits', [EEG.xmin EEG.xmax]*1000 , 'command', command, eegplotoptions{:}, 'butlabel', 'mark segments');
        % Now, "marked_segments" should contain start points and durations of manually marked segments (in pt)

        reply = input('After marking bad segments, press any key to continue: ','s');

        % Figure out good segments (in original data)
        % convert information of bad segments into information of good segments
        if ~isempty(marked_segments) 
            [good_seg_start, good_seg_dur, n_good_seg] = bad2good_segments(marked_segments.seg_start, marked_segments.seg_dur, EEG);

            % information of marked segments ready to be applied to original dataset
            good_segments.start = good_seg_start';
            good_segments.dur = good_seg_dur';
            EEG.etc.good_segments = good_segments; % put in EEG structure
            EEG_N20.etc.good_segments.start = EEG_N20.srate/down_sr * EEG.etc.good_segments.start; % put also in EEG structure for later anlyses; multiply by 20 since 250Hz * 20 = 5000Hz (srate)
            EEG_N20.etc.good_segments.dur = EEG_N20.srate/down_sr * EEG.etc.good_segments.dur;
        else % no segments were manually marked
            good_segments.start = 1;
            good_segments.dur = size(EEG.data, 2);
            EEG.etc.good_segments = good_segments; % put in EEG structure
            EEG_N20.etc.good_segments.start = 1; % set to 1 (all data)
            EEG_N20.etc.good_segments.dur = size(EEG_N20.data, 2); % all data 
        end 

    EEG_inspected = EEG; % copy inspected dataset
    
    
    if 0 %isfield(EEG.etc, 'channel_interpolation') % bypass visual data inspection
        EEG_marked = pop_loadset('filename',[subj '_sr250Hz_1to45Hz_vi_notremoved_OLD_without_channel_interpolation.set'],'filepath', savepath_full); % load old downsampled file with bad segment marks
        EEG.etc.good_segments = EEG_marked.etc.good_segments;
        EEG = pop_saveset(EEG, 'filename',[subj '_sr250Hz_1to45Hz_vi_notremoved.set'],'filepath', savepath_full); % save interpolated version (overwrite old one)
        
        EEG_N20.etc.good_segments.start = EEG_N20.srate/down_sr * EEG.etc.good_segments.start; % put also in EEG structure for later anlyses; multiply by 20 since 250Hz * 20 = 5000Hz (srate)
        EEG_N20.etc.good_segments.dur = EEG_N20.srate/down_sr * EEG.etc.good_segments.dur;
        
        EEG = pop_select(EEG, 'point', [EEG.etc.good_segments.start, EEG.etc.good_segments.start + EEG.etc.good_segments.dur]); % do not remove bad segments later, too
    
    % Usual Case (no bad channels):
    else 
        % save current preprocessing step before removing artifacts
        if ~exist(savepath_full)
            mkdir(savepath_full)
        end
        EEG = pop_saveset(EEG, 'filename',[subj '_sr250Hz_1to45Hz_vi_notremoved.set'],'filepath', savepath_full);

        % actually remove bad segments
        EEG_removed = eeg_eegrej(EEG,eegplot2event(TMPREJ, -1)); % 02_KA: 1440 events removed
        EEG = EEG_removed; % continue with corrected dataset
        clear EEG_removed EEG_loaded % do something good to the RAM
    end
    
    
    %% Re-reference to average reference
    % Get channel information
    EEG = pop_chanedit(EEG, 'lookup',electr_pos_file);
    % Insert FCz (original ref)
    EEG = pop_chanedit(EEG, 'insert', EEG.nbchan+1,'changefield',{EEG.nbchan+1 'labels' 'FCz'});
    EEG = pop_chanedit(EEG, 'setref',{['1:' num2str(EEG.nbchan+1)] 'FCz'});

    % EEG = pop_chanedit(EEG, 'insert',
    % EEG.nbchan+1,'changefield',{EEG.nbchan+1 'labels' 'POz'}); % for subject 02
    % EEG = pop_chanedit(EEG, 'setref',{['1:' num2str(EEG.nbchan+1)] 'POz'});  % for subject 02
    
    indxEog = find(ismember({EEG.chanlocs.labels}, {'VEOG1', 'VEOG2', 'HEOG1', 'HEOG2'})); % indices of EOG channels   
    EEG = pop_select(EEG,'nochannel',  indxEog);

    % Average reference
    EEG = pop_reref(EEG, [],'refloc',struct('labels',{'FCz'},'type',{''},'theta',{[]},'radius',{[]},'X',{[]},'Y',{[]},'Z',{[]},'sph_theta',{[]},'sph_phi',{[]},'sph_radius',{[]},'urchan',{[]},'ref',{'FCz'},'datachan',{0}));

    % % for subject 02:
    % EEG = pop_reref(EEG, [],'refloc',struct('labels',{'POz'},'type',{''},'theta',{[]},'radius',{[]},'X',{[]},'Y',{[]},'Z',{[]},'sph_theta',{[]},'sph_phi',{[]},'sph_radius',{[]},'urchan',{[]},'ref',{'POz'},'datachan',{0}));

    % save preproc dataset
    EEG = pop_saveset(EEG, 'filename',[subj '_sr250Hz_1to45Hz_vi_averef_notch.set'],'filepath', savepath_full);

    %% ICA
    % run ICA
    EEG = pop_runica(EEG, 'extended',1,'interupt','on'); % run ICA

    % save dataset containing ICA weights
    EEG = pop_saveset( EEG, 'filename',[subj '_sr250Hz_1to45Hz_vi_averef_notch_ICA.set'],'filepath', savepath_full);
    EEG_sr250Hz = EEG; % safe copy of dataset


    %% Preprocess dataset for N20 analyses (with sampling rate = 5kHz)
    % stimulation artifacts have already been interpolated; breaks etc have been removed; bad channels have been removed
    EEG = EEG_N20; % continue working with actual data for analysis
    eeglab redraw

    % remove bad segments (use prune parameters from above) -> above multiplied by EEG.srate/down_sr
    EEG = pop_select(EEG, 'point', [EEG.etc.good_segments.start, EEG.etc.good_segments.start + EEG.etc.good_segments.dur]); % take only experimental data
    EEG = pop_saveset(EEG, 'filename',[subj '_sr5kHz_before_filter.set'],'filepath', savepath_full);

    % band-pass filter
    filt_hz = [30 200];
    [b,a] = butter(2, filt_hz/(EEG.srate/2));
    EEG.data = filtfilt(b, a, double(EEG.data)')';

    % Re-referencing
    % Insert FCz (original ref)
    EEG = pop_chanedit(EEG, 'insert', EEG.nbchan+1,'changefield',{EEG.nbchan+1 'labels' 'FCz'});
    EEG = pop_chanedit(EEG, 'setref',{['1:' num2str(EEG.nbchan+1)] 'FCz'});
    % EEG = pop_chanedit(EEG, 'insert', EEG.nbchan+1,'changefield',{EEG.nbchan+1 'labels' 'POz'}); % subject 02
    % EEG = pop_chanedit(EEG, 'setref',{['1:' num2str(EEG.nbchan+1)] 'POz'}); % subject 02

    % Compute average reference (exclude EOG)
    indxEog = find(ismember({EEG.chanlocs.labels}, {'VEOG1', 'VEOG2', 'HEOG1', 'HEOG2'})); % indices of EOG channels
    EEG = pop_select(EEG,'nochannel',  indxEog);
    EEG = pop_reref(EEG, [], 'refloc',struct('labels',{'FCz'},'type',{''},'theta',{[]},'radius',{[]},'X',{[]},'Y',{[]},'Z',{[]},'sph_theta',{[]},'sph_phi',{[]},'sph_radius',{[]},'urchan',{[]},'ref',{'FCz'},'datachan',{0}));
    % EEG = pop_reref(EEG, [], 'refloc',struct('labels',{'POz'},'type',{''},'theta',{[]},'radius',{[]},'X',{[]},'Y',{[]},'Z',{[]},'sph_theta',{[]},'sph_phi',{[]},'sph_radius',{[]},'urchan',{[]},'ref',{'POz'},'datachan',{0})); % subject 02
    
    % Check channel location
    EEG = pop_chanedit(EEG, 'lookup',electr_pos_file);

    % Apply ICA decomposition matrices from preproc dataset to N20 dataset
    EEG = pop_editset(EEG, 'icachansind', 'EEG_sr250Hz.icachansind', 'chanlocs', '{EEG_sr250Hz.chanlocs EEG_sr250Hz.chaninfo EEG_sr250Hz.urchanlocs }', 'icaweights', 'EEG_sr250Hz.icaweights', 'icasphere', 'EEG_sr250Hz.icasphere');
    EEG = pop_saveset(EEG, 'filename',[subj '_sr5kHz_preproc_ICA.set'],'filepath', savepath_full);
    
    % Select bad components
    pop_selectcomps(EEG, [1:size(EEG.icaweights, 1)] );
    reply = input('After marking bad components, press any key to continue: ','s');
    
    % save components
    marked_comps_ICA = EEG.reject.gcompreject;
    save([savepath_full subj '_comp_marked_for_rejection.mat'], 'marked_comps_ICA') % save marked components

    % Remove bad ICs from N20 data
    EEG = pop_subcomp(EEG, [], 0);

    % Divide dataset into conditions    
    start_notask_lat = [EEG.event(  find(ismember({EEG.event.type}, num2str(start_notask_tr)), 1, 'last')  ).latency]; % get latency
    end_notask_lat = [EEG.event(  find(ismember({EEG.event.type}, num2str(end_notask_tr)), 1, 'last')  ).latency]; % get latency    
    start_task_lat = [EEG.event(  find(ismember({EEG.event.type}, num2str(start_task_tr)), 1, 'last')  ).latency]; % get latency
    end_task_lat = [EEG.event(  find(ismember({EEG.event.type}, num2str(end_task_tr)), 1, 'last')  ).latency]; % get latency
    
    EEG_notask = pop_select(EEG, 'point', [start_notask_lat - offset*EEG.srate, end_notask_lat + offset*EEG.srate]); % passive stimulation sequence
    EEG_task = pop_select(EEG, 'point', [start_task_lat - offset*EEG.srate, end_task_lat + offset*EEG.srate]); % somatosensory discrimination task
    
    % save dataset(s)
    EEG_notask = pop_saveset( EEG_notask, 'filename',[subj_names{s} '_notask_pchip_sr5kHz_' num2str(filt_hz(1)) 'to' num2str(filt_hz(2)) 'Hz_vi_averef_nonotch_ICA_removed.set'],'filepath', savepath_full);
    EEG_task = pop_saveset( EEG_task, 'filename',[subj_names{s} '_task_pchip_sr5kHz_' num2str(filt_hz(1)) 'to' num2str(filt_hz(2)) 'Hz_vi_averef_nonotch_ICA_removed.set'],'filepath', savepath_full);
    

    %% Preprocess peripheral data
    EEG = EEG_peri;    
    
    % Remove data segments as in EEG data
    EEG_eeg = pop_loadset('filename',[subj '_sr5kHz_before_filter.set'],'filepath', savepath_full);
    EEG.etc.good_segments = EEG_eeg.etc.good_segments;
    EEG = pop_select(EEG, 'point', [EEG.etc.good_segments.start, EEG.etc.good_segments.start + EEG.etc.good_segments.dur]); % take only experimental data
    
    % High-pass filter
    filt_hz = [70];
    [b,a] = butter(2, filt_hz/(EEG.srate/2), 'high');
    EEG.data = filtfilt(b, a, double(EEG.data)')'; % dataset for artifact detection and later ICA

    % Notch filters      
    [b,a] = butter(2, [48 52]/(EEG.srate/2), 'stop');
    EEG.data = filtfilt(b, a, double(EEG.data)')'; % dataset for artifact detection and later ICA
    
    [b,a] = butter(2, [148 152]/(EEG.srate/2), 'stop');
    EEG.data = filtfilt(b, a, double(EEG.data)')'; % dataset for artifact detection and later ICA
    
    % Divide dataset into conditions    
    start_notask_lat = [EEG.event(  find(ismember({EEG.event.type}, num2str(start_notask_tr)), 1, 'last')  ).latency]; % get latency
    end_notask_lat = [EEG.event(  find(ismember({EEG.event.type}, num2str(end_notask_tr)), 1, 'last')  ).latency]; % get latency    
    start_task_lat = [EEG.event(  find(ismember({EEG.event.type}, num2str(start_task_tr)), 1, 'last')  ).latency]; % get latency
    end_task_lat = [EEG.event(  find(ismember({EEG.event.type}, num2str(end_task_tr)), 1, 'last')  ).latency]; % get latency
    
    EEG_notask = pop_select(EEG, 'point', [start_notask_lat - offset*EEG.srate, end_notask_lat + offset*EEG.srate]); % passive stimulation sequence
    EEG_task = pop_select(EEG, 'point', [start_task_lat - offset*EEG.srate, end_task_lat + offset*EEG.srate]); % somatosensory discrimination task
    
    % save dataset(s)
    EEG_notask = pop_saveset( EEG_notask, 'filename',[subj_names{s} '_peri_notask_pchip_sr5kHz_HP' num2str(filt_hz(1)) 'Hz_vi_notch.set'],'filepath', savepath_full);
    EEG_task = pop_saveset( EEG_task, 'filename',[subj_names{s} '_peri_task_pchip_sr5kHz_HP' num2str(filt_hz(1)) 'Hz_vi_notch.set'],'filepath', savepath_full);  

end % subject loop














