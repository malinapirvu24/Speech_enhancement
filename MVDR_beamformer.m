function  S_hat = MVDR_beamformer (Rn, )

num_freq_bins = size(fft_frames, 1);
S_hat = zeros(num_freq_bins, num_frames);
S_hat_wiener = zeros(num_freq_bins, num_frames); % For Wiener filter estimate


for k = 1:num_freq_bins
    R_n_inv = pinv(R_n(:, :, k));  % Inverse of noise covariance matrix
    noise_power = mean(diag(R_W(:, :, k))); % Noise power estimation 
    for l = 1:num_frames
        Y = squeeze(fft_frames(k, l, :));  % Noisy signal for freq-bin k, frame l
        w_mvu = R_n_inv * ones(num_mics, 1) / (ones(1, num_mics) * R_n_inv * ones(num_mics, 1));  % MVU weights
        S_hat(k, l) = w_mvu' * Y;  % Clean signal estimate
    end
end

end 
