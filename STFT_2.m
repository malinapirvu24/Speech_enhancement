function STFT_signal = STFT_2(data, window_length, window_overlap,hop_size, impulse_response)

    num_mics = size(impulse_response,1);
 
    % Window creation 
    num_samples = length(data);
    num_windows = floor((num_samples-window_length)/hop_size)+1;
    windows = zeros(window_length, num_windows, num_mics);

    for m = 1:num_mics
        for n = 1:num_windows
            start_idx = (n-1)*hop_size + 1;
            end_idx = start_idx + window_length - 1;
            windows(:,n,m) = data(start_idx:end_idx).* impulse_response(m, :)';

        end
       
    end

    % FFT over the windows

    FFT_windows = fft(windows,[],1);

    % Keep positive frequencies and DC component
    STFT_signal = FFT_windows(1: ceil(window_length/2),:,:);



end