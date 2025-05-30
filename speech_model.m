% Clear environment
clear all
close all
clc

%% Import data

% Set the path to the folder containing the audio files
data_folder = '/Users/malinapirvu/Desktop/TuDelft/2024-2025/Q4/AP/Speech_enh_git';

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


%% STFT signals

% select window from: h_inter1, h_inter2, h_inter3, h_inter4, h_target
window = impulse_data.h_inter1;
window_length = size(window,2);
window_overlap = 0.5;
hop_size = window_length*(1-window_overlap);

% Import the source signals
s1 = audio_data.clean_speech_wav.signal;
s2 = audio_data.clean_speech_2_wav.signal;

% Perform convolution with a window ->FFT
% sizes of res_1 and res_2 are half of the window_length 
res_1 = STFT_2(s1,window_length,window_overlap,hop_size,window);
res_2 = STFT_2(s2,window_length,window_overlap,hop_size,window);

% Perform convolution with a window ->FFT for noise
noise = audio_data.babble_noise_wav.signal;
noise_STFT = STFT_2(noise,window_length,window_overlap,hop_size,window);

%% Padding to match dimensions 
% Get time dimensions
T_noise = size(noise_STFT, 2);
T_res1 = size(res_1, 2);
T_res2 = size(res_2, 2);

% Determine max time dimension
T_max = max([T_noise, T_res1, T_res2]);

% Pad each STFT matrix along time dimension (2nd dim)
% Padding with zeros at the end

noise_STFT_padded = padarray(noise_STFT, [0, T_max - T_noise, 0], 0, 'post');
res_1_padded      = padarray(res_1,      [0, T_max - T_res1,  0], 0, 'post');
res_2_padded      = padarray(res_2,      [0, T_max - T_res2,  0], 0, 'post');


X = res_1_padded + res_2_padded + noise_STFT_padded;


%% Compute noise covariance
num_noise_windows = floor((T_max - window_length) / hop_size) + 1;
num_mics = size(window,1);

R_n = zeros(num_mics, num_mics, size(noise_STFT_padded, 1));  % Frequency-bin specific covariance

for k = 1:size(noise_STFT_padded, 1)  % Loop over frequency bins
    for l = 1:num_noise_windows
        noise_vec = squeeze(noise_STFT_padded(k, l, :));  % Noise vector for freq-bin k, frame l
        R_n(:, :, k) = R_n(:, :, k) + (noise_vec * noise_vec') / num_noise_windows;
    end
end








