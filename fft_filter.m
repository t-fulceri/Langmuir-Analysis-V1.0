close all

cycle_select = 18;

V = smooth_voltage_cycles_rise(:,cycle_select);
I = smooth_current_cycles_rise(:,cycle_select);

%I curve
% plot(V,I,'.-');
% hold on;

%I' curve
scale1 = 1;
D_I = scale1*diff(I);
f_I_prime = figure;
movegui(f_I_prime,'northwest');
% plot(V(1:end-1),D_I,'.-');
plot(D_I,'.-');
title('dI/dV plot');

%Low-pass filtering through FFT

Fs = 10e6;                      % Sampling frequency (assume 1 MHz)
T = 1/Fs;                       % Sample time
L = length(I);                  % Length of signal
t = (0:L-1)*T;                  % Time vector

%Prepare the array of frequencies
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
% Use only positive frequencies as reference
freq = Fs/2*linspace(0,1,NFFT/2+1);

%Perform the Fast Fourier Transform using the array of frequencies
fft_D_I = fft(D_I,NFFT)/L;

% Plot only the low-frequency part of the spectrum
plot_end = round(length(freq)*1);

%For every +/- frequency pair, sum their contributions together
double_abs_fft_D_I = 2*abs(fft_D_I(1:NFFT/2+1));
%Find peak frequency in the plotted range
[M_fft_D_I,ind] = max(double_abs_fft_D_I(1:plot_end));

% Plot single-sided amplitude spectrum.
f_fft_D_I = figure;
movegui(f_fft_D_I,'north');
ax_fft_D_I = axes;
double_abs_fft_D_I_normalized = double_abs_fft_D_I./M_fft_D_I;
loglog(ax_fft_D_I,freq(1:plot_end),double_abs_fft_D_I_normalized(1:plot_end),'.-');
title(ax_fft_D_I,'Normalized dI/dV frequency spectrum');
xlabel(ax_fft_D_I,'Frequency [Hz]');
ylabel(ax_fft_D_I,'|I''(f)| = |FFT(I''(t))|');

%Extract Amplitude/Phase for each frequency
complex_fourier = zeros(length(freq),2);
complex_fourier(:,1) = abs(fft_D_I(1:NFFT/2+1));
complex_fourier(:,2) = angle(fft_D_I(1:NFFT/2+1));

time = linspace(1,L,L);
f_test = figure;
movegui(f_test,'east');
D_I_of_t = zeros(1,L);
for k = 1:length(complex_fourier);
    A_k = complex_fourier(k,1);
    phi_k = complex_fourier(k,2);
    D_I_of_t = D_I_of_t + A_k.*sin(f(k)*time + phi_k);
end

plot(time,D_I_of_t);

% %Kill high frequencies
% 
% %Estabilish cut-off freq
% bandwidth = length(fft_D_I)/2;
% cut_off = round(bandwidth/10);
% %Apply filter
% fft_D_I_aLP = fft_D_I;
% fft_D_I_aLP(cut_off:bandwidth) = 0;
% fft_D_I_aLP(bandwidth+cut_off-1:end) = 0;
% 
% %For every +/- frequency pair, sum their contributions together
% double_abs_fft_D_I_aLP = 2*abs(fft_D_I_aLP(1:NFFT/2+1));
% %Find peak frequency in the plotted range
% [M_fft_D_I_aLP,ind] = max(double_abs_fft_D_I_aLP(1:plot_end));
% 
% % Plot single-sided amplitude spectrum (after low-pass).
% f_fft_D_I_aLP = figure;
% movegui(f_fft_D_I_aLP,'center');
% ax_fft_D_I_aLP = axes;
% double_abs_fft_D_I_aLP_normalized = double_abs_fft_D_I_aLP./M_fft_D_I_aLP;
% loglog(ax_fft_D_I_aLP,freq(1:plot_end),double_abs_fft_D_I_aLP_normalized(1:plot_end),'.-');
% title(ax_fft_D_I_aLP,'Normalized dI/dV frequency spectrum after Low-Pass filter');
% xlabel(ax_fft_D_I_aLP,'Frequency [Hz]');
% ylabel(ax_fft_D_I_aLP,'|I''(f)| = |FFT(I''(t))|');