% Clear environment
clear all
close all
clc

%% Import data

% Set the path to the folder containing the audio files
%data_folder = '/Users/malinapirvu/Desktop/TuDelft/2024-2025/Q4/AP/Speech_enh_git';
data_folder = 'C:\workspace\matlab\tudelft\array_processing\speech_enhancement\Speech_enhancement';

% List of WAV files to read
wav_files = {
    'Speech_shaped_noise.wav', 
    'clean_speech.wav', 
    'clean_speech_2.wav', 
    'babble_noise.wav', 
    'aritificial_nonstat_noise.wav'
};

%% Read WAV files into a struct

% import 2 clean signals and 3 noise signals

audio_data = struct();
for k = 1:length(wav_files)
    file_path = fullfile(data_folder, wav_files{k});
    if exist(file_path, 'file')
        [y, Fs] = audioread(file_path);
        audio_data.(matlab.lang.makeValidName(wav_files{k})) = struct('signal', y, 'Fs', Fs);
    else
        warning('File not found: %s', file_path);
    end
end

%% Load impulse responses from .mat file
% the file contains 5 sets of impulse responses. 4 microphones and 1 source
% signal

impulse_path = fullfile(data_folder, 'impulse_responses.mat');
if exist(impulse_path, 'file')
    impulse_data = load(impulse_path);
else
    error('Impulse response file not found: %s', impulse_path);
end

%% Import all the sound signals
num_mics = 4; % multimicrohopne system with 4 mics

% For filtering get all impulse responses: h_inter1, h_inter2, h_inter3, h_inter4, h_target
h_inter1 = impulse_data.h_inter1;
h_inter2 = impulse_data.h_inter2;
h_inter3 = impulse_data.h_inter3;
h_inter4 = impulse_data.h_inter4;
h_target = impulse_data.h_target;

% Import the all the sound signals
s1 = audio_data.clean_speech_wav.signal;
s2 = audio_data.clean_speech_2_wav.signal;
noise_babble = audio_data.babble_noise_wav.signal;
noise_artificial = audio_data.aritificial_nonstat_noise_wav.signal;
noise_speech_shaped = audio_data.Speech_shaped_noise_wav.signal;

%% Filtering
% Select s1 as a target source signal (s2 and 3 noises are interferers + noise)
% Filter all the sound signals using the 5 impulse response sets
filter_length = size(h_target, 2); % all the impulse responses has same length

filtered_s1 = zeros(length(s1)+ filter_length - 1, num_mics);
filtered_s2 = zeros(length(s2)+ filter_length - 1, num_mics);
filtered_noise_babble = zeros(length(noise_babble)+ filter_length - 1, num_mics);
filtered_noise_artificial = zeros(length(noise_artificial)+ filter_length - 1, num_mics);
filtered_noise_speech_shaped = zeros(length(noise_speech_shaped)+ filter_length - 1, num_mics);

% Convolve all sound signal with corresponding impulse response
for m = 1:num_mics
    filtered_s1(:, m) = conv(s1, h_target(m, :));
    filtered_s2(:, m) = conv(s2, h_inter1(m, :));
    filtered_noise_babble(:, m) = conv(noise_babble, h_inter2(m, :));
    filtered_noise_artificial(:, m) = conv(noise_artificial, h_inter3(m, :));
    filtered_noise_speech_shaped(:, m) = conv(noise_speech_shaped, h_inter4(m, :));
end

%% STFT signals

% Window parameters
window_length = 0.02 *  Fs; % 20 ms frame length
window_overlap = 0.5; % 50% overlap
hop_size = window_length*(1-window_overlap);

% Perform convolution with a window ->FFT for clean signals
% sizes of res_1 and res_2 are half of the window_length 

res_s1 = STFT_2(filtered_s1, window_length, hop_size, num_mics); 
res_s2 = STFT_2(filtered_s2, window_length, hop_size, num_mics); 

% Perform convolution with a window ->FFT for noise
res_noise_babble = STFT_2(filtered_noise_babble, window_length, hop_size, num_mics); 
res_noise_artificial = STFT_2(filtered_noise_artificial, window_length, hop_size, num_mics);
res_noise_speech_shaped = STFT_2(filtered_noise_speech_shaped, window_length, hop_size, num_mics);

%% Padding to match dimensions 
%% Padding to match dimensions 

% Get num windows
num_windows_res1 = size(res_s1, 2);
num_windows_res2 = size(res_s2, 2);
num_windows_noise_babble = size(res_noise_babble, 2);
num_windows_noise_artificial = size(res_noise_artificial, 2);
num_windows_noise_speech_shaped = size(res_noise_speech_shaped, 2);

% Determine max num windows 
num_windows_padded = max([num_windows_res1, num_windows_res2, ...
    num_windows_noise_babble, num_windows_noise_artificial, ...
    num_windows_noise_speech_shaped]);

% Determine max num samples (padded)
num_samples_max = max([length(filtered_s1), length(filtered_s2), ...
    length(filtered_noise_babble), length(filtered_noise_artificial), ...
    length(filtered_noise_speech_shaped)]);

% Pad each STFT matrix along num windows (2nd dim)
% Padding with zeros at the end

res_s1_padded      = padarray(res_s1,      [0, num_windows_padded - num_windows_res1,  0], 0, 'post');
res_s2_padded      = padarray(res_s2,      [0, num_windows_padded - num_windows_res2,  0], 0, 'post');
res_noise_babble_padded = padarray(res_noise_babble, [0, num_windows_padded - num_windows_noise_babble, 0], 0, 'post');
res_noise_artificial_padded = padarray(res_noise_artificial, [0, num_windows_padded - num_windows_noise_artificial, 0], 0, 'post');
res_noise_speech_shaped_padded = padarray(res_noise_speech_shaped, [0, num_windows_padded - num_windows_noise_speech_shaped, 0], 0, 'post');

% Received signal
X = res_s1_padded + res_s2_padded + res_noise_babble_padded + ...
    res_noise_artificial_padded + res_noise_speech_shaped_padded;


%% Compute noise covariance
num_noise_windows = floor((num_samples_max - window_length) / hop_size) + 1;
num_freq_bins = size(X, 1);

R_n = zeros(num_mics, num_mics, num_freq_bins);  % Frequency-bin specific covariance

for k = 1:num_freq_bins  % Loop over frequency bins

    % Get interference (s2) and noise matrices
    s2_k = transpose(squeeze(res_s2_padded(k, :, :))); % Interference matrix for freq-bin k
    noise_babble_k = transpose(squeeze(res_noise_babble_padded(k, :, :)));  % Noise matrix for freq-bin k
    noise_artificial_k = transpose(squeeze(res_noise_artificial_padded(k, :, :))); 
    noise_speech_shaped_k = transpose(squeeze(res_noise_speech_shaped_padded(k, :, :))); 

    % Compute the individual covariances
    R_n_s2 = (s2_k * s2_k') / num_noise_windows;
    R_n_babble = (noise_babble_k * noise_babble_k') / num_noise_windows;
    R_n_artificial = (noise_artificial_k * noise_artificial_k') / num_noise_windows;
    R_n_speech_shaped = (noise_speech_shaped_k * noise_speech_shaped_k') / num_noise_windows;

    % Sum all the individual covariances
    R_n(:, :, k) = R_n(:, :, k) + (R_n_s2 + R_n_babble + R_n_artificial + R_n_speech_shaped) ;
end

%% Calculate acoustic transfer function and apply MVDR beamformer

% Calculate acoustic transfer function
A = calculate_ATF(X, R_n, num_mics, num_freq_bins);

res_hat = MVDR_beamformer(X, R_n, A, num_freq_bins, num_windows_padded);

%% Get estimated R_s and choose the reference microphone based on the average clean signal power across frequency bins 
R_s_hat = estimate_R_s(X, R_n, num_mics, num_freq_bins);

% Multi-channel Wiener beamformer
res_hat_Wiener = Wiener_beamformer(X, R_n, A, R_s_hat, num_freq_bins, num_windows_padded);

% Delay and sum beamformer
res_hat_DS = Delay_Sum_beamformer(X, A, num_freq_bins, num_windows_padded);
%% Reconstruct time-domain signal
S_recon_MVDR = zeros(num_samples_max,1);
S_recon_Wiener = zeros(num_samples_max,1);
S_recon_DS = zeros(num_samples_max,1);

% Does not need to divide by magnitude of the window cause its 1
%frame_weight_sum = zeros(num_samples, 1);  % Sum of Hann windows for overlap-add

for l = 1:num_windows_padded
    start_idx = (l - 1) * hop_size + 1;
    end_idx = start_idx + window_length-1 - 1;
    % For MVDR
    window_ifft = real(ifft([res_hat(:, l); conj(flip(res_hat(2:end, l), 1))]));  % IFFT to time-domain
    S_recon_MVDR(start_idx:end_idx) = S_recon_MVDR(start_idx:end_idx) + window_ifft; %does not need to divide by magnitude of the window cause its 1
    % For Wiener
    window_ifft_Wiener = real(ifft([res_hat_Wiener(:, l); conj(flip(res_hat_Wiener(2:end, l), 1))]));  % IFFT to time-domain
    S_recon_Wiener(start_idx:end_idx) = S_recon_Wiener(start_idx:end_idx) + window_ifft_Wiener; %does not need to divide by magnitude of the window cause its 1
    % For Delay and Sum
    window_ifft_DS = real(ifft([res_hat_DS(:, l); conj(flip(res_hat_DS(2:end, l), 1))]));  % IFFT to time-domain
    S_recon_DS(start_idx:end_idx) = S_recon_DS(start_idx:end_idx) + window_ifft_DS; %does not need to divide by magnitude of the window cause its 1
    %frame_weight_sum(start_idx:end_idx) = frame_weight_sum(start_idx:end_idx) + hann_window;
end
% Normalize by window sum
%S_recon = S_recon ./ frame_weight_sum;

% Save the result
%audiowrite('clean_signal_MVDR.wav', S_recon, Fs);

%% Beamformers performance evaluation using SNR_out

[SNR_MVDR, SNR_Wiener, SNR_DS] = SNR_out(R_n, A, R_s_hat, num_freq_bins);

mean_SNR_MVDR = mean(SNR_MVDR);
mean_SNR_Wiener = mean(SNR_Wiener);
mean_SNR_DS = mean(SNR_DS);

fprintf('Mean SNR (MVDR): %.2f dB\n', 10*log10(mean_SNR_MVDR));
fprintf('Mean SNR (Wiener): %.2f dB\n', 10*log10(mean_SNR_Wiener));
fprintf('Mean SNR (DS): %.2f dB\n', 10*log10(mean_SNR_DS));

%% Beamformers performance evaluation using MSE

[MSE_MVDR, MSE_Wiener, MSE_DS] = calculate_MSE(R_n, A, R_s_hat, num_freq_bins);

fprintf('Mean MSE (MVDR): %.2f \n', mean(MSE_MVDR));
fprintf('Mean MSE (Wiener): %.2f \n', mean(MSE_Wiener));
fprintf('Mean MSE (DS): %.2f \n', mean(MSE_DS));







