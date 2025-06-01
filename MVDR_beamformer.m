function  res_hat = MVDR_beamformer(X, R_n, A, num_freq_bins, num_windows_padded)

res_hat = zeros(num_freq_bins, num_windows_padded);


for k = 1:num_freq_bins
    R_n_inv = pinv(R_n(:, :, k));  % Inverse of noise covariance matrix
    for l = 1:num_windows_padded
        x = squeeze(X(k, l, :));  % Noisy signal for freq-bin k, frame l
        a = transpose(A(k, :)); %vector acoustic transfer func for freq bin k, dth source signal
        w_mvdr = R_n_inv * a / (a' * R_n_inv * a);  % MVU weights
        res_hat(k, l) = w_mvdr' * x;  % Clean dth source signal estimate 
    end
end






