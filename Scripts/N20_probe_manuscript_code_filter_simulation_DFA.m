%%% Temporal signatures of critical dynamics in cortical excitability as probed by early somatosensory responses %%%
% Authors: Tilman Stephani, Gunnar Waterstraat, Stefan Haufe, Garbiel Curio, Arno Villringer & Vadim V. Nikulin (2020)

%%% Code for simulations on influence of temporal filtering on DFA exponents %%%
% -> look for power-law dynamics in filtered stochastically independent data


%% Create pseudo SEP (varying stochastically across trials)
EEG = pop_loadset('/data/pt_01972/Preproc_data/N20_study1/12_TS/12_TS_notask_pchip_sr5kHz_30to200Hz_vi_averef_nonotch_ICA_removed.set'); % take shape of the SEP from empirical data
EEG = pop_epoch(EEG, {'A - Out'}, [-0.1 0.6]);
n_trials = EEG.trials;
n_pnts = EEG.pnts;

mean_SEP = mean(EEG.data(42,:,:),3);
mean_SEP = (mean_SEP-mean(mean_SEP)) / mean(mean_SEP) * (-1); % normalize: mean=0; var=1

amp_factor = rand(n_trials, 1); % random amplification of every trial
sig = mean_SEP' * amp_factor'; % now, SEP should vary randomly
    
    % plot average
    %figure, plot(EEG.times, mean(sig,2))
    
    % plot single trials
    %figure, imagesc(sig')
    
sig = reshape(sig, n_pnts*n_trials, 1);
    

%% generate pink noise with same length and filter it
offset = 20; % try different start points for filtering
noise = pinknoise(length(sig)+offset); % mean=0, var=1
noise = noise(offset+1:end);

filt_hz = [30 200];
[b,a] = butter(2, filt_hz/(sr/2));
noise = filtfilt(b, a, double(noise));

%     [P,f]= pwelch(noise,hanning(sr),0,sr,sr);
%     figure; semilogy(f,P)
%     title('pink noise')
    
% mix signal and noise
SNR = 2;
sig_mix = SNR*sig + noise';
sig_mix = reshape(sig_mix, n_pnts, n_trials);

% trial-wise DFA
start = 7;
stop = 70;
num_segment=30;
start_fit = 7;
stop_fit  = 70;
fg=0; 
prune = [1 size(sig_mix,2)-1]; 

exponent = [];
parfor i = 1:size(sig_mix,1)
    [exponent(i),Amplitude,Alpha,time,st,epochs]=dfa_2018(sig_mix(i,:), start,stop,num_segment,start_fit,stop_fit,prune, fg);
end
   
% plot DFA time course
figure, plot(EEG.times, exponent)





    
    