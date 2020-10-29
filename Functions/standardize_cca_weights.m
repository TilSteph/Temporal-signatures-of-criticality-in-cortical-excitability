function [cca_weights_standardized, iCCA_tangential, iCCA_radial] = standardize_cca_weights(CCA_comps, cca_weights, Ast, EEG_orig_epoched)
% Standardize the first two CCA component weights (polarity)
% 1) Choose polarity of CCA weights in such a way, that the N20 is a negative peak
% 2) "Classify" as tangential or radial -> does not work yet!
% Author: Tilman Stephani, Jan 2019

%%
% %% Inputs (for debugging)
% EEG_orig = pop_loadset('filename',[subj '_sr5kHz_30to200Hz_vi_averef_NEWnotch_ICA_removed.set'],'filepath', [datapath_pre subj_names{s} '/']);
% EEG_orig.data = double(EEG_orig.data);
% 
% [CCA_comps, cca_weights, Ast, EEG_orig_epoched, ix_accepted] = apply_CCA(EEG_orig, {'A - Out', 'B - Out'}, [5 80], [-100 200], 2);


%%
% %% Check signals
% %  Spatial patterns of CCA components
% figure;
% subplot(1,2,1); topoplot(Ast(:,1), EEG_orig_epoched.chanlocs); title('CCA 1')
% subplot(1,2,2); topoplot(Ast(:,2), EEG_orig_epoched.chanlocs); title('CCA 2')
% 
% % Averages of CCA components
% x_lim = [-5 100];
% figure;
% subplot(1,2,1)
% plot(EEG_orig_epoched.times, mean(CCA_comps(1,:,:), 3))
% xlabel('time (ms)')
% ylabel('component activation (a.u.)')
% xlim(x_lim)
% title('CCA component 1')
% 
% subplot(1,2,2)
% plot(EEG_orig_epoched.times, mean(CCA_comps(2,:,:), 3))
% xlabel('time (ms)')
% ylabel('component activation (a.u.)')
% xlim(x_lim)
% title('CCA component 2')
% 
% % Original SEP at posterior sites
% elec = find(ismember({EEG_orig_epoched.chanlocs.labels}, 'P4'));
% figure;
% plot(EEG_orig_epoched.times, mean(EEG_orig_epoched.data(elec, :, :),3))
% xlim(x_lim)
% title(EEG_orig_epoched.chanlocs(elec).labels)


%% Correlate sensor space signal with CCA components (averages)
%elec = find(ismember({EEG_orig_epoched.chanlocs.labels}, 'P4')); % choose electrode in sensor space
elec = find(ismember({EEG_orig_epoched.chanlocs.labels}, 'CP4')); % choose electrode in sensor space


% averages
Sensor_average = squeeze(mean(EEG_orig_epoched.data(elec, :, :),3));
CCA1_average = squeeze(mean(CCA_comps(1,:,:), 3));
CCA2_average = squeeze(mean(CCA_comps(2,:,:), 3));
CCA3_average = squeeze(mean(CCA_comps(3,:,:), 3));

% correlations
[rho1, pval1] = corr(Sensor_average', CCA1_average');
[rho2, pval2] = corr(Sensor_average', CCA2_average');
[rho3, pval3] = corr(Sensor_average', CCA3_average');

%% Flip CCA weight polarity in case of negative correlation
fprintf(['\nChoosing polarity for CCA component 1 & 2 of subject ' EEG_orig_epoched.subject ' ... \n' ...
         'rho(' EEG_orig_epoched.chanlocs(elec).labels ',CCA1) = ' num2str(rho1) '; p = ' num2str(pval1) '\n' ...
         'rho(' EEG_orig_epoched.chanlocs(elec).labels ',CCA2) = ' num2str(rho2) '; p = ' num2str(pval2) '\n' ...
         'rho(' EEG_orig_epoched.chanlocs(elec).labels ',CCA3) = ' num2str(rho3) '; p = ' num2str(pval3) '\n']);

cca_weights_standardized = cca_weights;
if rho1 < 0
    %CCA_comps(1,:,:) = CCA_comps(1,:,:) * (-1); % flip polarity of component
    cca_weights_standardized(:,1) = cca_weights(:,1) * (-1);
    fprintf('Polarity of CCA1 has been flipped. \n');
else
    fprintf('Keep original polarities of CCA1. \n');
end

if rho2 < 0
    %CCA_comps(2,:,:) = CCA_comps(2,:,:) * (-1); % flip polarity of component
    cca_weights_standardized(:,2) = cca_weights(:,2) * (-1);
    fprintf('Polarity of CCA2 has been flipped. \n');
else
    fprintf('Keep original polarities of CCA2. \n');
end

if rho3 < 0
    %CCA_comps(2,:,:) = CCA_comps(2,:,:) * (-1); % flip polarity of component
    cca_weights_standardized(:,3) = cca_weights(:,3) * (-1);
    fprintf('Polarity of CCA3 has been flipped. \n\n');
else
    fprintf('Keep original polarities of CCA3. \n\n');
end

% Warnings if correlations are low
if pval1 > 0.05, fprintf(['   Warning: low correlation for CCA1 (p=' num2str(pval1) ') \n']); end
if pval2 > 0.05, fprintf(['   Warning: low correlation for CCA2 (p=' num2str(pval2) ') \n']); end
if pval3 > 0.05, fprintf(['   Warning: low correlation for CCA3 (p=' num2str(pval3) ') \n\n']); end



%% Classify as tangential or radial components (first two components)
elec_front = find(ismember({EEG_orig_epoched.chanlocs.labels}, 'Fz'));
elec_back = find(ismember({EEG_orig_epoched.chanlocs.labels}, 'P4'));
elecs_surround = find(ismember({EEG_orig_epoched.chanlocs.labels}, {'F4', 'Fz', 'Pz', 'PO4', 'FT8', 'TP8'}));
elec_centr = find(ismember({EEG_orig_epoched.chanlocs.labels}, 'CP4'));

Dip_diff = Ast(elec_front,1:4) - Ast(elec_back,1:4); % calculate difference between back and front electrode
[~, iCCA_tangential] = max(abs(Dip_diff)); % search for "dipole" component

Rad_diff = Ast(elec_centr,1:4) - mean(Ast(elecs_surround,1:4));
[~, iCCA_radial] = max(abs(Rad_diff)); % search for radial component

end

