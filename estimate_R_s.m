function R_s_hat = estimate_R_s(X, R_n, num_mics, num_freq_bins)
%single target source model

R_s_hat = zeros(num_mics, num_mics, num_freq_bins);

for k = 1:num_freq_bins
    % Get noise covariance and data matrix X corresponding to particular freq. k
    R_n_k = R_n(:, :, k); % size (num_mics) x (num_mics)
    X_k = squeeze(X(k, :, :)); % size n x num_mics
    % Since R_n corresponding to each freq is symmetric, we can apply EVD
    [V, D] = eig(R_n_k);
    
    % If different noises are added to each microphone
    R_n_inv_half = V * diag(1./sqrt(diag(D))) * V'; %R_n^-1/2
    R_n_half = V .* sqrt(D) * V; %R_n^1/2

    % Pre-whiten data
    X_whiten = R_n_inv_half * transpose(X_k); % size (num_mics) x n
    % Compute covariance of X_whiten
    R_x_whiten = (X_whiten * X_whiten') / size(X_whiten, 2);
    % Apply EVD
    [U, D] = eig(R_x_whiten);
    [eigvals_sorted, idx] = sort(diag(D), 'descend');
    
    U_modified = U(:, idx(1));
    
    % Estimate R_s
    lambda = eigvals_sorted(1);
    R_s_whiten_hat = lambda * (U_modified * U_modified'); 
    R_s_hat(:, :, k) = real(R_n_half * R_s_whiten_hat * R_n_half);

end
