%%% Temporal signatures of critical dynamics in cortical excitability as probed by early somatosensory responses %%%
% Authors: Tilman Stephani, Gunnar Waterstraat, Stefan Haufe, Garbiel Curio, Arno Villringer & Vadim V. Nikulin (2020)

%%% Code for generating individual head models in Brainstorm %%%

% Based on a pipeline by Stefan Haufe; adapted by Tilman Stephani and Alice Hodapp (2019)

% General processing steps:
%        1. Import MRI
%        2. Automatically adjust fiducials (needs to be visually checked by
%        figures). These are stored in the "subjectimage_T1.mat" and the
%        information about the fiducials: "MRI_back.SCS.NAS" and
%        "MRI_back.NCS.AC"
%        3. Import the freesurfer folder and overwrite the Fiducials
%        4. Compute the BEM (skull thickness might be adjusted)
%        5. Import electrodes coordinates file:
%            - a "fake" eeg file needs to be imported
%            - import electrode coordinate file (not aligned)
%            - project electrodes on scalp (creates individualized coordinate file for each subject)
%        6. Compute Head Model



%% environment preparation
clear 

% brainstorm software needs to be open to run the functions of brainstorm:
addpath('/data/p_01972/Brainstorm/brainstorm3/');
brainstorm


%% GUI: make new protocol
protocol_name = 'N20_final_corrected';

%% Define names and paths
db_path = ['/data/pt_01972/Brainstorm/brainstorm_db/' protocol_name '/'];
mri_path = '/data/pt_01972/Coregistration/Freesurfer_output/';
funeeg_path = '/data/pt_01972/EEG_means/';
digitize_path = '/data/pt_01972/Coregistration/Renamed_digitize_files/';

subj_dir = dir('/data/p_01972/RAW_EEG/N20_study1/*_*');
subj_names = {subj_dir.name};
subj_in = [2:27, 29:33]; % exclude subjects without individual MRI scans (for these, the head models were created via the Brainstorm GUI)
all_ID = {subj_names{subj_in}}; 
N=length(all_ID);

for i=1:N
subj = all_ID{i};
    
digitize_file = dir([digitize_path 'channel_' subj '*']);
digitize_file = [digitize_file.folder '/' digitize_file.name];

%% Define T1 files

sFiles = [];
SubjectNames = {...
    subj};

RawFiles = {...
    [mri_path subj '/mri/T1.mgz']}; 

% Start a new report
bst_report('Start', sFiles);

% Process: Import MRI
sFiles = bst_process('CallProcess', 'process_import_mri', sFiles, [], ...
    'subjectname', SubjectNames{1}, ...
    'mrifile',     {RawFiles{1}, 'MGH'}, ...
    'nas',         [0, 0, 0], ...
    'lpa',         [0, 0, 0], ...
    'rpa',         [0, 0, 0], ...
    'ac',          [0, 0, 0], ...
    'pc',          [0, 0, 0], ...
    'ih',          [0, 0, 0]);

% Process: Compute MNI transformation (SET Fiducials)
sFiles = bst_process('CallProcess', 'process_mni_affine', sFiles, [], ...
    'subjectname', SubjectNames{1});


% we load the information about the Fiducials which are stored in subjects folder here:
MRI_back = load([db_path '/anat/' subj '/subjectimage_T1.mat']);

%% LOAD ANATOMY (FREESURFER) FOLDER  
RawFiles = {...
    [mri_path subj], ...
    [funeeg_path 'EEG_mean_' subj '.set'], ... % it's not necessary to link raw file, so directly import EEG Average 
    digitize_file};

% Process: Import anatomy folder and overwrite the fiducials 
sFiles = bst_process('CallProcess', 'process_import_anatomy', sFiles, [], ...
    'subjectname', SubjectNames{1}, ...
    'mrifile',     {RawFiles{1}, 'FreeSurfer'}, ...
    'nvertices',   5000, ... % urspr?nglich auf 15000, 5000 reichen aber
    'nas',         MRI_back.SCS.NAS, ...
    'lpa',         MRI_back.SCS.LPA, ...
    'rpa',         MRI_back.SCS.RPA, ...
    'ac',          MRI_back.NCS.AC, ...
    'pc',          MRI_back.NCS.PC, ...
    'ih',          MRI_back.NCS.IH);
 

%% Set mid cortex as default surface 

%get file
SurfaceFile = [subj '/tess_cortex_mid_low.mat']
[sSubject, iSubject, iSurface] = bst_get('SurfaceFile', SurfaceFile);

%change default 
db_surface_default(iSubject, 'Cortex', iSurface);

% Save database
db_save(1); %1=always save protocol, too

%% Process: Generate BEM surfaces
sFiles = bst_process('CallProcess', 'process_generate_bem', sFiles, [], ...
    'subjectname', SubjectNames{1}, ...
    'nscalp',      1922, ...
    'nouter',      1922, ...
    'ninner',      1922, ...
    'thickness',   4);


%% EEG and channel file
% Process: Import EEG avarage 
sFiles = bst_process('CallProcess', 'process_import_data_time', sFiles, [], ...
    'subjectname',  SubjectNames{1}, ...
    'condition',    '', ...
    'datafile',     {{RawFiles{2}}, 'EEG-EEGLAB'}, ...
    'timewindow',   [-0.1, 0.08], ...
    'split',        0, ...
    'ignoreshort',  0, ...
    'channelalign', 0, ...
    'usectfcomp',   1, ...
    'usessp',       1, ...
    'freq',         [], ...
    'baseline',     []);

% Process: Set channel file (should contain the electrodes ordered in the same way as the EEG dataset!)
sFiles = bst_process('CallProcess', 'process_import_channel', sFiles, [], ...
    'channelfile',  {RawFiles{3}, 'BST'}, ...
    'usedefault',   1, ...  % NotAligned: GSN HydroCel 128 E1
    'channelalign', 0, ...
    'fixunits',     0, ... % if 1, electrodes flipped
    'vox2ras',      1);

% Process: Refine registration
for q = 1:3 % make sure that best fit is found
    bst_process('CallProcess', 'process_headpoints_refine', sFiles, []);
end
% Process: Project electrodes on scalp
bst_process('CallProcess', 'process_channel_project', sFiles, []);

end

%% COMPUTE HEAD MODEL
for i=1:N
subj = all_ID{i};

sFiles = [];
sFiles = [subj '/EEG_mean_' subj '/data_block001.mat'];
FileType = file_gettype(sFiles);

% Process: Compute head model
sFiles = bst_process('CallProcess', 'process_headmodel', sFiles, [], ...
    'Comment',     '', ...
    'sourcespace', 1, ...  % Cortex surface
    'volumegrid',  struct(...
         'Method',        'isotropic', ...
         'nLayers',       17, ...
         'Reduction',     3, ...
         'nVerticesInit', 4000, ...
         'Resolution',    0.005, ...
         'FileName',      ''), ...
    'meg',         1, ...  % 
    'eeg',         3, ...  % OpenMEEG BEM
    'ecog',        1, ...  % 
    'seeg',        1, ...  % 
    'openmeeg',    struct(...
         'BemFiles',     {{}}, ...
         'BemNames',     {{'Scalp', 'Skull', 'Brain'}}, ...
         'BemCond',      [1, 0.0125, 1], ...
         'BemSelect',    [1, 1, 1], ...
         'isAdjoint',    1, ... % falls es abst?rzt aufgrund von Memory Problemen hier auf 1 setzen
         'isAdaptative', 1, ...
         'isSplit',      1, ...
         'SplitLength',  4000));

end 







