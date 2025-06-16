function  [MSE_MVDR, MSE_Wiener, MSE_DS] = calculate_MSE(R_n, A, R_s1, num_freq_bins, num_mics)

MSE_MVDR = zeros(num_freq_bins, 1);
MSE_Wiener = zeros(num_freq_bins, 1);
MSE_DS = zeros(num_freq_bins, 1);

for k = 1:num_freq_bins
    R_n_inv = pinv(R_n(:, :, k));  % Inverse of noise covariance matrix
    R_n_k = R_n(:, :, k);
    R_s_k = R_s1(:, :, k);

    %Choose reference mic as the one with highest average power
    [~, ref_mic_index] = max(diag(R_s_k));
    var_s_ref = R_s_k(ref_mic_index, ref_mic_index);
    corr_s_ref = zeros(num_mics, 1); % cross-covariance vector between s (in all mics) and refernce s
    corr_s_ref(ref_mic_index) = var_s_ref; 

    a = transpose(A(k, :)); %vector acoustic transfer func for freq bin k, dth source signal

    % Compute beamformers
    %w_single_channle_Wiener = var_s_ref / ( var_s_ref + 1/((a' * R_n_inv * a)));
    w_single_channle_Wiener = 0.0029 / ( 0.0029 + 1/((a' * R_n_inv * a)));
    w_mvdr = R_n_inv * a / (a' * R_n_inv * a);  % MVU weights
    w_wiener = w_single_channle_Wiener * w_mvdr; % Wiener weights
    w_ds =  a / (a' * a);  % Delay and Sum weights

    % Calculate SNR_out for all beamformers
    MSE_MVDR(k) = real((w_mvdr' * R_s_k * w_mvdr) - (w_mvdr' * corr_s_ref) - (corr_s_ref' * w_mvdr) + (w_mvdr' * R_n_k * w_mvdr));
    MSE_Wiener(k) = real((w_wiener' * R_s_k * w_wiener) - (w_wiener' * corr_s_ref) - (corr_s_ref' * w_wiener) + (w_wiener' * R_n_k * w_wiener));
    MSE_DS(k) = real((w_ds' * R_s_k * w_ds) - (w_ds' * corr_s_ref) - (corr_s_ref' * w_ds) + (w_ds' * R_n_k * w_ds));
end






