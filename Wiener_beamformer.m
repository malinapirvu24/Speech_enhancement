function  res_hat = Wiener_beamformer(X, R_n, A, R_s_hat, num_freq_bins, num_windows_padded)

res_hat = zeros(num_freq_bins, num_windows_padded);


for k = 1:num_freq_bins
    R_n_inv = pinv(R_n(:, :, k));  % Inverse of noise covariance matrix
    R_s_k = R_s_hat(:, :, k);
    %Choose reference mic as the one with highest average power
    [~, ref_mic_index] = max(diag(R_s_k));
    var_s_ref = R_s_k(ref_mic_index, ref_mic_index);
    for l = 1:num_windows_padded
        x = squeeze(X(k, l, :));  % Noisy signal for freq-bin k, frame l
        a = transpose(A(k, :)); %vector acoustic transfer func for freq bin k, dth source signal
        w_single_channle_Wiener = var_s_ref / ( var_s_ref + 1/((a' * R_n_inv * a)));
        w_mvdr = R_n_inv * a / (a' * R_n_inv * a);  % MVU weights
        w = w_single_channle_Wiener * w_mvdr;
        res_hat(k, l) = w' * x;  % Clean dth source signal estimate 
    end
end






