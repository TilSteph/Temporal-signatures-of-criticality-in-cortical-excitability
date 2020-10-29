function [CCA_comps, Wst, Ast, R, iCCA_tangential, iCCA_radial, EEG_epoched, ix_accepted] = apply_CCA(EEG, epoch_markers, train_win, eval_win, num_comp, flag_standardize)
% Apply canonical correlation analysis to EEG data (eeglab format).
%   Inputs:
%   EEG: EEGLAB dataset with EEG.data in format 'channel by time' (continuous EEG data)
%   epoch_markers: markers which are used for epoching
%   train_win: window for which CCA weights are found; in ms
%   eval_win: window to which CCA weights are applied; in ms
%   num_comp: number of components to be extracted
%   flag_standardize: option for polarity standardization
%   
%   Outputs:
%   CCA_comps: CCA components (component by time by epoch)
%   Wst: CCA weights
%   Ast: spatial patterns for all components (pattern by component number)
%   R: canonical correlation of components with original data
%   EEG_epoched: epoched input dataset
%   ix_accepted: indices of epochs that were long enough for epoching
%
%   EEGLAB should be initialized before running this function.
% Author: Tilman Stephani, 11/2018

% % debug
% epoch_markers = {'A - Out', 'B - Out'};
% eval_win = [-100 600];
% train_win = [5 80];
% num_comp = 4;

fprintf('\nDo CCA ... \n')

%% Epoch data
epoch_length_t = eval_win/1000;
[EEG, ix_accepted] = pop_epoch(EEG, epoch_markers, epoch_length_t, 'epochinfo', 'yes'); % use second output later for choosing epochs (e.g. in behavioral data) AND prune parameters in DFA!!!!
    %EEG = pop_rmbase(EEG, [-30 0]); % in case of low HP filter!!!
    
%% Prepare long and short versions of the dataset
EEG_long = EEG;
EEG_short = pop_select(EEG,'time', train_win/1000);

%% Prepare concatenated matrices
% average matrix
ERP_mean = mean(EEG_short.data(:, :, :),3)';        
Mav_short = repmat(ERP_mean, EEG_short.trials, 1);

% single trials
Mst_short = reshape(EEG_short.data, EEG_short.nbchan, EEG_short.pnts*EEG_short.trials)';
Mst_long = reshape(EEG_long.data, EEG_long.nbchan, EEG_long.pnts*EEG_short.trials)';
    
%% Do CCA
[Wav, Wst, R]=canoncorr(Mav_short, Mst_short);
%n_cca_all = size(R_all, 2);
n_cca_keep = num_comp;

% Apply obtained weights to long dataset (according to eval_win)
%Wst(:,2) = Wst(:,2) * (-1); % only test!!
CCA_concat = (Mst_long * Wst(:,1:num_comp))';

% Spatial patterns
Ast = cov(Mst_short)*Wst;

%% Re-reshape
CCA_comps = reshape(CCA_concat, num_comp, EEG.pnts, EEG.trials);


%% output epoched EEG
EEG_epoched = EEG_long;

%% Polarity standardization & component classification (tangential vs radial)
if flag_standardize % if chosen
    [Wst, iCCA_tangential, iCCA_radial] = standardize_cca_weights(CCA_comps, Wst, Ast, EEG_epoched); % standardize polarity in such a way that the N20 in CCA space is a negative peak
    
    % Re-calculate
    CCA_concat = (Mst_long * Wst(:,1:num_comp))';
    CCA_comps = reshape(CCA_concat, num_comp, EEG.pnts, EEG.trials);
    Ast = cov(Mst_short)*Wst;
end

end

