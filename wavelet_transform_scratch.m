function [data_wavelet_freq, data_fft, wavelet_fft, frex] = wavelet_transform_scratch(time_series_data, fs, frex, cycles, wavelet_half_win_size)  
% Wavelet transform code based on Mike X Cohen text book chapter 13.11
% Yunshu Fan 2020-05-11
% inputs: 
%  time_series_data = ft_data.trial{1}; % can be any time-series data. 
%       NOTE: rows and columns cannot be swapped: Rows: data of each channel; Columns: data of each time point
%       If only one channel, this can be a 1 x N vector 
%  fs = data.fs % sampling frequency of time_series_data
%  freq_space_type = 'log'; % 'log' or 'lin' for log- or linear spacing
%  min_freq =  1; % lower bound freq
%  max_freq = 40; % upper bound freq
%  num_frex = 20; % number of freqs
%  num_cycle_type = 'var'; % 'var' or 'fix'
%  num_cycle = [3, 10]; % 1 element if num_cycle_type = 'var'; 2 elements if num_cycle_type = 'var'
%  wavelet_half_win_size = 1; % in unit of seconds. decide how long the
%  wavelet is. The wavelet with the lowest frequency should taper to ~0
%  within this window. If not sure, try with the function plot_wavelet.m
%
% outputs:
%   data_wavelet_freq: wavelet transform of the time series data. 
%       dimensions 1~3: number of frequencies x number of time points x number of channels
%   data_fft: frequency domain of data. 
%       dimensions 1, 2: number of frequencies  x n_conv_pow2 (n_conv_pow2
%       is for the convenience of fft
%   wavelet_fft: frequency domain of wavelet. 
%       dimensions 1, 2: number of frequencies  x n_conv_pow2 (n_conv_pow2
%       is for the convenience of fft
%   frex: a list of frequencies of the wavelets

%  check if the length of the frequency and cycle vectors matches
if length(frex)~=length(cycles)
    error('length and frequency and cycle vectors must match')
end

% define wavelet parameters
time = -wavelet_half_win_size:1/fs:wavelet_half_win_size;
s = cycles./(2*pi*frex);

% definte convolution parameters
n_wavelet            = length(time);
n_data               = size(time_series_data, 2);
n_convolution        = n_wavelet+n_data-1;
n_conv_pow2          = pow2(nextpow2(n_convolution));
half_of_wavelet_size = (n_wavelet-1)/2;

% loop through frequencies and compute wavelet
num_frex = length(frex);
wavelet_fft = zeros(num_frex, n_conv_pow2);
for fi=1:num_frex
    % no normalization term for gaussian
    wavelet_time_domain = exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2)));
%     % with normalization term (wrong term, what's used in the text book)
%     wavelet_time_domain = sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2)));
    cmwX = fft(wavelet_time_domain, n_conv_pow2);
    cmwX = cmwX / max(cmwX);
    wavelet_fft(fi, :) = cmwX;
end

% get FFT of data
% data_fft = (fft(time_series_data', n_conv_pow2))'; % fft is performed to each column, but I want the return to be each channel in each row
% data_fft = flip((fft(time_series_data', n_conv_pow2))'); % fft is performed to each column, but I want the return to be each channel in each row
data_fft = fft(time_series_data, n_conv_pow2, 2); % fft is performed to each column, but I want the return to be each channel in each row

% initialize
[nchannels, ntimepoints] = size(time_series_data);
data_wavelet_freq = zeros(num_frex, ntimepoints, nchannels);
% loop through channels and freq
for fi = 1:num_frex
    for ichannel = 1:nchannels
        % convolution
        eegconv = ifft(wavelet_fft(fi, :).*data_fft(ichannel, :));
        eegconv = eegconv(1:n_convolution);
        data_wavelet_freq(fi, :, ichannel) = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size);
        fprintf('finished channel %d/%d freq %d/%d\n', ichannel, nchannels, fi, num_frex);
    end
end