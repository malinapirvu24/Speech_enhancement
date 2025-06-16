function  res_hat = Delay_Sum_beamformer(X, A, num_freq_bins, num_windows_padded)

res_hat = zeros(num_freq_bins, num_windows_padded);


for k = 1:num_freq_bins
    for l = 1:num_windows_padded
        x = squeeze(X(k, l, :));  % Noisy signal for freq-bin k, frame l
        a = transpose(A(k, :)); %vector acoustic transfer func for freq bin k, dth source signal
        w_ds =  a / (a' * a);  % Delay and Sum weights
        res_hat(k, l) = w_ds' * x;  % Clean dth source signal estimate 
    end
end






