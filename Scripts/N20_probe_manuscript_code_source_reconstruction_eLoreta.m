%%% Temporal signatures of critical dynamics in cortical excitability as probed by early somatosensory responses %%%
% Authors: Tilman Stephani, Gunnar Waterstraat, Stefan Haufe, Garbiel Curio, Arno Villringer & Vadim V. Nikulin (2020)

%%% Code for applying eLoreta to reconstruct EEG sources %%%


%% environment preparation
addpath('/data/pt_01972/Brainstorm/brainstorm3/');
brainstorm

addpath('/data/p_01972/EEGLAB/eeglab14_1_1b/');
eeglab
pop_editoptions( 'option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 1, 'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 0, 'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0, 'option_checkversion', 1, 'option_chat', 0);

addpath('/data/p_01972/Scripts/Code_data_manuscript_N20_power_law_dynamics/Functions/')

path_toolbox = '/data/p_01972/Other_toolboxes/';
addpath(genpath([path_toolbox,'tensor_toolbox_2.6']));
%addpath('/data/p_01972/Scripts/eLoreta')

% Note: in the osf repository, '/Functions/' contains all functions needed
% for eLoreta (apart from the tensor toolbox (commercial))


%% Define paths and names
%protocol_name = 'N20_final_all_subjects';
protocol_name = 'N20_final_corrected'; % use from now on! (but DFA_S computed on N20_final_all_subjects)
db_path = '/data/pt_01972/Brainstorm_Tilman/Brainstorm_db_Tilman/';

data_path = [db_path protocol_name '/data/'];
anat_path = [db_path protocol_name  '/anat/'];

savepath_pre = '/data/pt_01972/Preproc_data/N20_study1/'; % first part of save directory
datapath_pre = '/data/pt_01972/Preproc_data/N20_study1/'; % first part of data directory
savepath_figures = '/data/pt_01972/Results/figures/';

subj_dir = dir([anat_path '/*_*']); % here, brainstorm database is used to get subject IDs!
subj_names = {subj_dir.name};
% subj_in = [1:12, 14:17, 19:33]; % exclude subject 13
% subj_names = {subj_names{subj_in}}; 

% define conditions
cond = 'notask'; % passive stimulation sequence

% load CCA data 
A_all_tangential = [];
for s = 1:33
    
    datapath_full = [datapath_pre subj_names{s} '/']; % complete save directory     
    pl_suffx = ['_' cond '_30to200Hz_nonotch_stdzd_CP4'];
    
    CCA = load([datapath_full 'CCA_' subj_names{s} pl_suffx '.mat']); % CCA data including patterns
    load([datapath_pre 'Classify_CCA_task_30to200Hz_nonotch_stdzd_CP4.mat']); % load classification of components
    
    A_all_tangential(:,s) = CCA.Ast(:, comps_tangential(s)); % spatial pattern of the tangential CCA component
    disp(num2str(s))
end



%% subject loop
OutputFiles = {};
for s = 1:33     
    disp(subj_names{s})    

    %% Load individual headmodel
    hm = load([data_path subj_names{s} '/@default_study/headmodel_surf_openmeeg.mat']); 

    % Constrain sources (only dipoles normal to the cortex)
    Gain_constr = bst_gain_orient(hm.Gain, hm.GridOrient);
    
    %% eLORETA
    % apply eLoreta
    N_Ch = size(A_all_tangential,1); % # of channels
    
    H = eye(N_Ch) - ones(N_Ch) ./ N_Ch; % this operator centralizes the matrix which is mltiplied to it
    Gain_constr = double(ttm(tensor(Gain_constr),H,1)); % centralized leadfiled
    gamma = 0.05; % eLORETA reg. factor, This value was suggested by Stefan and Vadim.
    A_eloreta = mkfilt_eloreta_v2(Gain_constr,gamma); % eLORETA demixing matrix

    %% Reconstruct sources of sensor space data
    X_sensor = A_all_tangential(:,s); % tangential CCA
    S = A_eloreta' * X_sensor;
    
     % make sure default surface is "mid..."
    if ~ismember(subj_names{s}, {'01_BT', '03_PA', '27_TM', '31_MG'}) %s~=1 && s~=3 && s~=27 && s~=31 % individual brains
        SurfaceFile = [subj_names{s} '/tess_cortex_mid_low.mat'];
        [sSubject, iSubject, iSurface] = bst_get('SurfaceFile', SurfaceFile);
        db_surface_default(iSubject, 'Cortex', iSurface); %change default 
        db_save(1); %1=always save protocol, too    
    elseif ismember(subj_names{s}, {'27_TM'}) % s==27
        SurfaceFile = [subj_names{s} '/tess_cortex_pial_low.mat'];
        [sSubject, iSubject, iSurface] = bst_get('SurfaceFile', SurfaceFile);
        db_surface_default(iSubject, 'Cortex', iSurface); %change default 
        db_save(1); %1=always save protocol, too
    end
    
    %% create brainstorm results file
    results_struct = struct;
    results_struct.DataFile =  [];
    results_struct.Comment = 'Sources_CCA_tangential_eLoreta_task';
    %results_struct.Comment = 'Sources_CCA_radial_eLoreta';
    results_struct.HeadModelFile = [data_path subj_names{s} '/@default_study/headmodel_surf_openmeeg.mat'];
    results_struct.HeadModelType = 'surface';
    results_struct.History = [];
    results_struct.ImageGridAmp = S;
    results_struct.ImagingKernel = [];
    results_struct.Options = [];
    results_struct.SurfaceFile = [subj_names{s} '/tess_cortex_mid_low.mat'];
    results_struct.Time = 1;
    results_struct.nAvg = 1;
    results_struct.nComponents = 1;
    
    if ismember(subj_names{s}, {'01_BT', '31_MG'}) %s==1 || s==31 % standard anatomies
        results_struct.SurfaceFile = '@default_subject/tess_cortex_pial_high_5000V.mat';
    elseif ismember(subj_names{s}, {'03_PA'}) %s==3
        results_struct.SurfaceFile = '03_PA/tess_cortex_pial_low_5000V.mat';
    elseif ismember(subj_names{s}, {'27_TM'}) %s==27 % "cortex" surface not "mid"
        results_struct.SurfaceFile = '27_TM/tess_cortex_pial_low.mat';    
    end

    % save in database
    [sStudy, iStudy, iData] = bst_get('DataFile', ['/' subj_names{s} '/EEG_mean_NEW_' subj_names{s} '/data_block001.mat']); % get iStudy
    if isempty(sStudy), [sStudy, iStudy, iData] = bst_get('DataFile', ['/' subj_names{s} '/EEG_mean_' subj_names{s} '/data_block001.mat']); end % get iStudy
    results_struct_file = db_add(iStudy, results_struct); % save and register in database

    %% project to template surface    
    destSurfFile = '@default_subject/tess_cortex_pial_high_90000V.mat';
    OutputFiles{s} = bst_project_sources({results_struct_file}, destSurfFile, 0, 0); 
        
end
save([savepath_pre 'CCA_tangential_projected_on_template_eLoreta_filenames_final_corrected.mat'], 'OutputFiles')


%% Average patterns
load_ix = load([savepath_pre 'CCA_tangential_projected_on_template_eLoreta_filenames_final_corrected.mat']);
OutputFiles = load_ix.OutputFiles;

subj_in = [1:12, 14:17, 19:33]; % exclude subjects 13 and 18
sFiles = {}; ix_files = 0;
for ss = subj_in %1:length(OutputFiles) % do some anoying formatting
    ix_files = ix_files +1;
    sFiles{ix_files} = cell2mat(OutputFiles{ss});
end

% Process: Average: Everything (Brainstorm)
sFiles = bst_process('CallProcess', 'process_average', sFiles, [], ...
    'avgtype',         1, ...  % Everything
    'avg_func',        2, ...  % mean(abs); before: 1 arithmetic average: mean(x)
    'weighted',        0, ...
    'scalenormalized', 0);


% --> Reconstructed sources were visualized using the Brainstorm GUI.





