%%% Temporal signatures of critical dynamics in cortical excitability as probed by early somatosensory responses %%%
% Authors: Tilman Stephani, Gunnar Waterstraat, Stefan Haufe, Garbiel Curio, Arno Villringer & Vadim V. Nikulin (2020)

%%% Code for SNR simulations %%%


%% prepare environment
clear
addpath('/data/p_01972/Scripts/Code_data_manuscript_N20_power_law_dynamics/Functions/');
savepath_MANUSCRIPT = '/data/p_01972/Scripts/Code_data_manuscript_N20_power_law_dynamics/Data/'; % save path for preprocessed data that can be published according to the European GDPR


%% simulation parameters
true_DFA_list = 0.5 : 0.005 : 0.8; % vary DFA exponents of the signal
ntrials = 1000; % number of simulated trials (according to our empirical data)
len_signal = 1500;
n_distr = 100; % simulate every parameter combination several times to account for randomness in signal generation


%% create signals
rng('shuffle')
test_sig = linspace(0.001,6,len_signal); % continuously increasing signal; this the SNR is continuously manipulated
noise  = repmat(rand(1, ntrials, n_distr), len_signal,1,1); % create white noise (same noise vector for all SNR levels)


%% DFA parameters (according to the DFA applied to our empirical data)
start = 7;
stop = 70;
num_segment = 20;
start_fit = 7;
stop_fit  = 70;
fg=0;
prune = [1 ntrials-1];

exponent_noise = zeros(n_distr, 1);
exponent_true_emp = zeros(length(true_DFA_list), n_distr);
exponent_mix = zeros(len_signal, length(true_DFA_list), n_distr);
rms_noise = zeros(n_distr, 1);
SNR_test = zeros(len_signal, n_distr);
SNR_test_powerlaw = zeros(len_signal, length(true_DFA_list), n_distr);

%% Calculate SNR and DFA exponents of noise (should be close to 0.5)
parfor j = 1:n_distr
    % Calculate noise energy
    rms_noise(j) = rms(noise(1,:,j),2); % rms of white noise
    SNR_test(:, j) = test_sig / rms_noise(j);
    
    % DFA noise
    [exponent_noise(j),Amplitude,Alpha,time,st,epochs]=dfa_2018(squeeze(noise(1,:,j)), start,stop,num_segment,start_fit,stop_fit,prune, fg);
end


%% Mix signal and noise for all combinations of DFA exponents and SNR
parfor k = 1:length(true_DFA_list)
    rng(4); % arbitrary seed for rand()
    ERP_modulation = simulate_powerlaw(ntrials, true_DFA_list(k)); % generate a time course that exhibits a certain power-law relationship
    test_sig_powerlaw = test_sig' * (ERP_modulation+1)'; % convolve test signal with power-law time course; add 1 to modulation parameter (to make it positive) -> make whole signal obey a power-law across trials

    for n = 1:n_distr % repeat mixing procedure several times (to account for randomness in noise generation)
        % Calculate SNR (in case power-law conversion changes it)
        SNR_test_powerlaw(:,k,n) = mean(test_sig_powerlaw,2) / rms_noise(n);

        % Mix test signal with white noise (at all SNR levels)
        test_sig_mix = test_sig_powerlaw + noise(:,:,n);
        
        % Calculate DFA exponents of test signal and of mixed signal
        [exponent_true_emp(k,n),Amplitude,Alpha,time,st,epochs]=dfa_2018(ERP_modulation, start,stop,num_segment,start_fit,stop_fit,prune, fg);

        for m = 1:len_signal % DFA exponents for every SNR level
            [exponent_mix(m,k,n),Amplitude,Alpha,time,st,epochs]=dfa_2018(squeeze(test_sig_mix(m,:)), start,stop,num_segment,start_fit,stop_fit,prune, fg);
        end
    end    
    disp(k)
end

% save results
save([savepath_MANUSCRIPT 'DFA_mixing_distribution.mat'], 'exponent_mix', 'exponent_noise', 'exponent_true_emp', 'true_DFA_list', 'test_sig', 'rms_noise', 'SNR_test', 'SNR_test_powerlaw', 'noise')



