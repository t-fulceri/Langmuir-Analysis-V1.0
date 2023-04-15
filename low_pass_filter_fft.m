function [ filtered_signal ] = low_pass_filter_fft( signal, cutoff_fraction, freq_sampling )
%LOW PASS FILTER FFT Applies a low pass filter with spceified cut-off frequency
%   Tiziano Fulceri 2018-01-31
    
    %Fraction must be <= 1
    if cutoff_fraction > 1
        return
    end
    % Sampling frequency
    switch nargin
        case 2
            Fs = 1e7;  %DEFAULT VALUE (10 MHz) 
        case 3
            Fs = freq_sampling;
    end
    
    % Length of signal
    L = max(size(signal));             
    %Prepare the array of frequencies
    NFFT = 2^nextpow2(L); % Next power of 2 from length of y
    % Use only positive frequencies as reference
    freq = Fs/2*linspace(0,1,NFFT/2+1);
    %Perform the Fast Fourier Transform using the array of frequencies
    fft_signal = fft(signal,NFFT)/L;
    %Consider only positive frequencies
    fft_positive = fft_signal(1:NFFT/2+1);
    %For every +/- frequency pair, sum their amplitudes
    double_abs_fft_signal = 2*abs(fft_positive);
    
    
    %Prepare output frequencies array
    freq_cutoff = cutoff_fraction*freq_sampling;
    %The calculated freq_cutoff must be approximated by a frequency contained in the freq array:
    [min_distance, cutoff_index] = min(abs(freq - freq_cutoff));
    %Prepare output frequencies array
    freq_array = freq(1:cutoff_index);
    %Preapre output phases array
    phase_array = angle(fft_positive(1:cutoff_index));
    
    %Reconstruction of the filtered signal
    %Have to adjust by removing the mean value...
    %Not really sure why it is necessary, but it works
    filtered_signal = -mean(signal);
    time_step = (1/Fs);
    t = time_step*linspace(0,L-1,L);
    for n = 1:cutoff_index
        filtered_signal = filtered_signal + double_abs_fft_signal(n)*cos(2*pi*freq_array(n)*t + phase_array(n));
    end
    
    
%     %PLOT (JUST FOR DEBUG)
%     % Plot original signal plus extracted main component
%     f_supimp = figure;
%     movegui(f_supimp,'center');
%     ax_supimp = axes;
%     hold on;
%     plot(signal);
%     plot(filtered_signal);
%     hold off;
%     title(ax_supimp,'Original signal superimposed with low-frequency components only');
%     xlabel(ax_supimp,'Time [timesteps]');
%     ylabel(ax_supimp,'Phase');

end

