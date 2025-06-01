function STFT_signal = STFT_2(data, window_length, hop_size, num_mics)
    
    % Window creation 
    num_samples = length(data);
    num_windows = floor((num_samples-window_length)/hop_size)+1;
    windows = zeros(window_length, num_windows, num_mics);
    hann_window = hann(window_length, 'periodic');
    
    for m = 1:num_mics
        for n = 1:num_windows
            start_idx = (n-1)*hop_size + 1;
            end_idx = start_idx + window_length - 1;
            windows(:,n,m) = data(start_idx:end_idx, m).* hann_window; 
        end       
    end

    % FFT over the windows

    FFT_windows = fft(windows,[],1);

    % Keep positive frequencies and DC component
    STFT_signal = FFT_windows(1: ceil(window_length/2),:,:);

end