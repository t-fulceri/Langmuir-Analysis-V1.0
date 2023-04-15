close all

V = smooth_voltage;
I = smooth_current;

Fs = 10e6;                      % Sampling frequency (assume 1 MHz)
T = 1/Fs;                       % Sample time
L = length(smooth_voltage);     % Length of signal
t = (0:L-1)*T;                  % Time vector

%Voltage
f_V = figure;
movegui(f_V,'northwest');
ax_V = axes;
plot(ax_V,Fs*t,V,'.-b');
title(ax_V,'Bias voltage over time');
xlabel(ax_V,'time [arbitrary units]');
ylabel(ax_V,'V [V]');

%Current
f_I = figure;
movegui(f_I,'west');
ax_I = axes;
plot(ax_I,Fs*t,I,'.-r');
title(ax_I,'Probe current over time');
xlabel(ax_I,'time [arbitrary units]');
ylabel(ax_I,'I [A]');


%Prepare the array of frequencies
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
% Use only positive frequencies as reference
freq = Fs/2*linspace(0,1,NFFT/2+1);

%Perform the Fast Fourier Transform using the array of frequencies
fft_V = fft(V,NFFT)/L;

% Plot only the low-frequency part of the spectrum
plot_end = round(length(freq)*1);

%For every +/- frequency pair, sum their contributions together
double_abs_fft_V = 2*abs(fft_V(1:NFFT/2+1));
%Find peak frequency in the plotted range
[M_fft_V,ind] = max(double_abs_fft_V(1:plot_end));

% Plot single-sided amplitude spectrum.
f_fftV = figure;
movegui(f_fftV,'north');
ax_fftV = axes;
double_abs_fft_V_normalized = double_abs_fft_V./M_fft_V;
loglog(ax_fftV,freq(1:plot_end),double_abs_fft_V_normalized(1:plot_end),'.-b');
title(ax_fftV,'Normalized bias voltage frequency spectrum');
xlabel(ax_fftV,'Frequency [Hz]');
ylabel(ax_fftV,'|V(f)| = |FFT(V(t))|');

%Perform the Fast Fourier Transform using the array of frequencies
fft_I = fft(I,NFFT)/L;

%For every +/- frequency pair, sum their contributions together
double_abs_fft_I = 2*abs(fft_I(1:NFFT/2+1));
%Find peak frequency in the plotted range
[M_fft_I,ind] = max(double_abs_fft_I(1:plot_end));

% Plot single-sided amplitude spectrum.
f_fftI = figure;
movegui(f_fftI,'center');
ax_fftI = axes;
double_abs_fft_I_normalized = double_abs_fft_I./M_fft_I;
loglog(ax_fftI,freq(1:plot_end),double_abs_fft_I_normalized(1:plot_end),'.-r');
title(ax_fftI,'Normalized probe current frequency spectrum');
xlabel(ax_fftI,'Frequency [Hz]');
ylabel(ax_fftI,'|I(f)| = |FFT(I(t))|');

%delta_amplitude vector:
% double_abs_delta = double_abs_fft_I_normalized - double_abs_fft_V_normalized;
% f_delta_amplitude = figure;
% movegui(f_delta_amplitude,'south');
% ax_delta_amplitude = axes;
% semilogx(ax_delta_amplitude,freq(1:plot_end),double_abs_delta(1:plot_end),'.-k');
% title(ax_delta_amplitude,'Normalized amplitude difference between current and voltage');
% xlabel(ax_delta_amplitude,'Frequency [Hz]');
% ylabel(ax_delta_amplitude,'\DeltaA');


% Plot V phases.
V_phase = angle(fft_V(1:NFFT/2+1));
f_phaseV = figure;
movegui(f_phaseV,'northeast');
ax_phaseV = axes;
semilogx(ax_phaseV,freq(1:plot_end),V_phase(1:plot_end),'.-b');
ax_phaseV.YTick = [ -pi 0 pi ];
ax_phaseV.YTickLabel = ({ '-\pi','0','\pi' });
title(ax_phaseV,'Bias voltage phase');
xlabel(ax_phaseV,'Frequency [Hz]');
ylabel(ax_phaseV,'\phi_V');

% Plot I phases.
I_phase = angle(fft_I(1:NFFT/2+1));
f_phaseI = figure;
movegui(f_phaseI,'east');
ax_phaseI = axes;
semilogx(ax_phaseI,freq(1:plot_end),I_phase(1:plot_end),'.-r');
ax_phaseI.YTick = [ -pi 0 pi ];
ax_phaseI.YTickLabel = ({ '-\pi','0','\pi' });
title(ax_phaseI,'Probe current phase');
xlabel(ax_phaseI,'Frequency [Hz]');
ylabel(ax_phaseI,'\phi_I');


%Generate Delta_phi vector:
delta_phi = I_phase - V_phase;
f_delta_phi = figure;
movegui(f_delta_phi,'southeast');
ax_delta_phi = axes;
semilogx(ax_delta_phi,freq(1:plot_end),delta_phi(1:plot_end),'.-k');
ax_delta_phi.YTick = [ -pi 0 pi ];
ax_delta_phi.YTickLabel = ({ '-\pi','0','\pi' });
title(ax_delta_phi,'Phase difference between current and voltage');
xlabel(ax_delta_phi,'Frequency [Hz]');
ylabel(ax_delta_phi,'\Delta\phi = (\phi_I - \phi_V)');

%Find peak frequency in the plotted range
[M_fft_V,ind] = max(double_abs_fft_V(1:plot_end));
peak_freq = freq(ind)
delta_phi_at_freq_peak = delta_phi(ind)

%Idea: Calculate the fft of impedance:
fft_Z = fft_V./fft_I;
double_abs_fft_Z = 2*abs(fft_Z(1:NFFT/2+1));
f_fftZ = figure;
movegui(f_fftZ,'south');
ax_fftZ = axes;
loglog(ax_fftZ,freq(1:plot_end),double_abs_fft_Z(1:plot_end),'.-k');
title(ax_fftZ,'Impedance frequency spectrum');
xlabel(ax_fftZ,'Frequency [Hz]');
ylabel(ax_fftZ,'|I(f)/V(f)| = |Z(f)| = |FFT(Z(t))|');