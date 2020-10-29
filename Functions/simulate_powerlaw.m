function X = simulate_powerlaw(nr_samples, alpha)
% creates time series specified alpha-exponent, algorithm according to:
% Kasdin, N. Jeremy.
% "Discrete simulation of colored noise and stochastic processes
% and 1/f^\alpha power law noise generation."
% Proceedings of the IEEE 83.5 (1995): 802-827.

% IN:
%     nr_samples [integer] : length of time series
%     alpha      [integer] : time series has specified alpha-exponent
% OUT:
%       X [nr_samples x 1] : time series with specified alpha-exponent

beta = 2*alpha-1;
Q_d = 1;          %white noise will be in the range [-Q_d,Q_d]

%  generate the coefficients h_k.
hfa = zeros(2*nr_samples, 1);
hfa(1) = 1;
for i = 2:nr_samples
    hfa(i) = hfa(i-1) * (beta/2+(i - 2))/(i-1);
end

% fill the sequence w_k with white noise and pad with zeroes
wfa = [-Q_d + 2 * Q_d.* rand(nr_samples, 1); zeros(nr_samples, 1)];

% perform the discrete Fourier transforms
fh = fft(hfa);
fw = fft(wfa);

% multiply the two complex vectors and pad with zeroes
complex_prod = fh .* fw;
complex_prod = [complex_prod(1:nr_samples+1); zeros(nr_samples-1,1)];

%  inverse Fourier transform the result.
X = ifft(complex_prod);
X = real(X(1:nr_samples));

end