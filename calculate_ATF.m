function A = calculate_ATF(X, R_n, num_mics, num_freq_bins)
%single target source model

A = zeros(num_freq_bins, num_mics);
%R_s_hat = zeros(num_mics, num_mics);

for k = 1:num_freq_bins
    % Get noise covariance and data matrix X corresponding to particular freq. k
    R_n_k = R_n(:, :, k); % size (num_mics) x (num_mics)
    X_k = squeeze(X(k, :, :)); % size n x n
    % Since R_n corresponding to each freq is symmetric, we can apply EVD
    [V, D] = eig(R_n_k);

    % [eigenvals, idx] = sort(diag(D), 'descend');
    % V1 = V(:, idx(1));
    % D1 = eigenvals(1); % take largest eigenvalue (others are 0)
    % 
    % R_n_inv_half = 1/ sqrt(D1) * (V1 * V1'); %R_n^-1/2
    % R_n_half = sqrt(D1) * (V1 *  V1'); %R_n^1/2
    
    % If different noises are added to each microphone
    R_n_inv_half = V * 1./sqrt(diag(D)) * V'; %R_n^-1/2
    R_n_half = V .* sqrt(diag(D)) * V; %R_n^1/2

    % Pre-whiten data
    X_whiten = R_n_inv_half * transpose(X_k); % size (num_mics) x n
    % Compute covariance of X_whiten
    R_x_whiten = (X_whiten * X_whiten') / size(X_whiten, 2);
    % Apply EVD
    [U, D] = eig(R_x_whiten);
    [~, idx] = sort(diag(D), 'descend');
    
    U_modified = U(:, idx(1));
    
    % % Estimate R_s
    % lambda = largest_eigenvals(d);
    % R_s_whiten_hat = lambda * U_modified(:, d) * U_modified(:, d)'; 
    % R_s_hat(:, :, d) = R_n-half * R_s_whiten_hat * R_n_half;
    
    % Compute acoustic transfer function for single target
    A(k, :) = R_n_half * U_modified;
end
