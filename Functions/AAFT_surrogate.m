function surrogate_set = AAFT_surrogate(X, nreps)
% generate a surrogate dataset with same spectrum and distribution as X
% adapted: Tilman Stephani, 08/2019

% Code taken from: N Schaworonkow, DAJ Blythe, J Kegeles, G Curio, VV Nikulin: 
% Power-law dynamics in neuronal and behavioral data introduce spurious 
% correlations. Human Brain Mapping. 2015.
% http://doi.org/10.1002/hbm.22816

X = X(:);
nr_samples = numel(X);

surrogate_set = zeros(nr_samples, nreps);
for i = 1:nreps
    
    if mod(i,1000)==0
        display(['run iteration: ' num2str(i) '/' num2str(nreps)])
    end
    
    % create white noise vector with n entries
    white_noise = sort(randn(nr_samples,1));
    % Sort z and extract the ranks
    [sorted_z, ranks] = sort(X);
    [~, idx_ranks] = sort(ranks);
    
    % random phase surrogate on white noise
    y = fft(white_noise(idx_ranks)');
    Z_amps = abs(y);
    rand_phases = rand(1,floor(nr_samples/2))*2*pi;
    start = length(rand_phases);
    
    if mod(nr_samples,2) == 0
        start = start-1;
    end
    
    rand_phases = [0, rand_phases, -rand_phases(start:-1:1)];
    % put amps and phases together for complex Fourier spectrum
    white_noise_rand_phase = Z_amps .* exp(1i*rand_phases);
    % project the complex spectrum back to the time domain
    white_noise_rand_phase = real(ifft(white_noise_rand_phase));
    
    % extract the ranks of the phase randomized white noise
    [~, ranks] = sort(white_noise_rand_phase);
    [~, idx_ranks] = sort(ranks);
    
    % assign ranks of phase randomized normal deviates to sorted_z,
    % obtain AAFT surrogates
    z_surrogate = sorted_z(idx_ranks);
    z_surrogate = z_surrogate(:);
    
    surrogate_set(:,i) = z_surrogate;
    
end

end