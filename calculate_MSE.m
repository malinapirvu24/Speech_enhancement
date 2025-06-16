function  [MSE_MVDR, MSE_Wiener, MSE_DS] = calculate_MSE(R_n, A, R_s1, num_freq_bins)

MSE_MVDR = zeros(num_freq_bins, 1);
MSE_Wiener = zeros(num_freq_bins, 1);
MSE_DS = zeros(num_freq_bins, 1);

for k = 1:num_freq_bins
    R_n_inv = pinv(R_n(:, :, k));  % Inverse of noise covariance matrix
    R_n_k = R_n(:, :, k);
    R_s_k = R_s1(:, :, k);

    %Choose reference mic as the one with highest average power
    [~, ref_mic_index] = max(diag(R_s_k));
    
    corr_s_ref = R_s_k(:,ref_mic_index);
    var_s_ref = R_s_k(ref_mic_index, ref_mic_index);

    a = transpose(A(k, :)); %vector acoustic transfer func for freq bin k, dth source signal

    % Compute beamformers
    w_single_channle_Wiener = var_s_ref / ( var_s_ref + 1/((a' * R_n_inv * a)));
    w_mvdr = R_n_inv * a / (a' * R_n_inv * a);  % MVU weights
    w_wiener = w_single_channle_Wiener * w_mvdr; % Wiener weights
    w_ds =  a / (a' * a);  % Delay and Sum weights

    % Calculate SNR_out for all beamformers
    MSE_MVDR(k) = real(w_mvdr' * (R_s_k + R_n_k) * w_mvdr + var_s_ref - 2 * real(w_mvdr' * corr_s_ref));
    MSE_Wiener(k) = real(w_wiener' * (R_s_k + R_n_k) * w_wiener + var_s_ref - 2 * real(w_wiener' * corr_s_ref));
    MSE_DS(k) = real(w_ds' * (R_s_k + R_n_k) * w_ds + var_s_ref - 2 * real(w_ds' * corr_s_ref));
end






