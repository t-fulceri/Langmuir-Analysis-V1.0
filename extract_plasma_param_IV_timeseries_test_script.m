%Clear workspace
clear all;
%Close all figures
close all;

%Directory and filepath:

% %ILYA STUFF
% meas_directory = '\\pc-e5-ws-2\D\Ilya\04112016\500\';
% filename = 'meas_0770.h5';
% voltage_sign = 1;
% current_sign = 1;

% % %TIZIANO STUFF
% meas_directory = '\\pcvineta3\Data\tzf\LP\'
% filename = 'meas_0243.h5'
% voltage_sign = +1;
% current_sign = -1; %(IV-Curve convention: invert current axis);
% %Probe tip area
% L_tip = 8.00*1e-3; %8.00 mm length
% r_tip = 0.25*1e-3; %0.25 mm radius
% %Result:
% S_tip = L_tip*pi*r_tip;
% % %Ion mass number
% A = 39.948;%Argon

% %TIZIANO STUFF
meas_directory = '\\pc-e5-ws-2\D\tzf\Plasma_Profile_050kHz_4V_5_rep_CentralZone_REPRISE\'
filename = 'meas_0303.h5';
voltage_sign = +1;
current_sign = -1; %(IV-Curve convention: invert current axis);
%Probe tip area
L_tip = 4.00*1e-3; %4.00 mm length
r_tip = 0.25*1e-3; %0.25 mm radius
S_tip = L_tip*pi*r_tip;
% Ion mass number
A = 39.948;%Argon

%SEARCH PLASMA GUN ELECTRODE VOLTAGE RECORDING-----------------------------
%Add directories for functions
addpath('\\pcvineta3\Data\MATLAB\files');
addpath('\\pcvineta3\Data\MATLAB\files\HDF5-tools');

%Reconnection monitor output files directory
recmon_output_dir = '\\pc-e5-ws-2\D\reconnection_monitor\';

recmon_output_filename = 'rcc_219319.h5';

recmon_output_full_file_path = strcat(recmon_output_dir,recmon_output_filename);

[out,timing]=readRecMon(recmon_output_full_file_path, 'rogCalib.txt');
%--------------------------------------------------------------------------


[ return_status, output_extract_plasma_param_IV_timeseries ] = extract_plasma_param_IV_timeseries( meas_directory, filename, voltage_sign, current_sign, S_tip, A)
if(return_status ~= 0)
    fprintf('extract_plasma_param_IV_timeseries failed to extract plasma parameters\n');
end

%Floating potential
V_float = output_extract_plasma_param_IV_timeseries.V_float;
%Ion saturation current
I_sat_i_1st = output_extract_plasma_param_IV_timeseries.I_sat_i_1st;
%Temperature
T_e_lin_fit = output_extract_plasma_param_IV_timeseries.T_e_lin_fit;

T_e_min = output_extract_plasma_param_IV_timeseries.T_e_min;
T_e_max = output_extract_plasma_param_IV_timeseries.T_e_max;
T_e_median = output_extract_plasma_param_IV_timeseries.T_e_median;

%Electron density
n_e_derived = output_extract_plasma_param_IV_timeseries.n_e_derived;

n_e_min = output_extract_plasma_param_IV_timeseries.n_e_min;
n_e_max = output_extract_plasma_param_IV_timeseries.n_e_max;
n_e_median = output_extract_plasma_param_IV_timeseries.n_e_median;

aesthetic_range = [250 700];

time_axis = output_extract_plasma_param_IV_timeseries.time_axis;
%Move along the time axis untill the beginning of the plasma shot
time_axis =  time_axis - timing.recMon.start*1e6;

total_time_axis = output_extract_plasma_param_IV_timeseries.total_time_axis;
%Move along the time axis untill the beginning of the plasma shot
total_time_axis = total_time_axis - timing.recMon.start*1e6;

%Plot
f_v = figure;
%Avoid graphical errors
set(gcf,'Renderer','painters');

time_label_string = 'Time [µs]';

movegui(f_v,'northwest');
subplot(3,1,1);
plot((out(1).t)*1e6, out(1).data);
title('Electron Gun 1 Cathode Curent');
xlabel(time_label_string);
ylabel('I_{Cathode} [A]');

%time_limits = xlim;
time_limits = aesthetic_range;
xlim(time_limits);




subplot(3,1,2);
plot(total_time_axis,output_extract_plasma_param_IV_timeseries.voltage,'-');
title('Voltage over time');
xlabel(time_label_string);
ylabel('V [V]');
xlim(time_limits);
subplot(3,1,3);
plot(total_time_axis,output_extract_plasma_param_IV_timeseries.current,'-');
title('Current over time');
xlabel(time_label_string);
ylabel('I [A]');
xlim(time_limits);



%Plot the superimposed voltage cycles
f = figure;
set(gcf,'Renderer','painters');

movegui(f,'north');
subplot(2,1,1);
plot(output_extract_plasma_param_IV_timeseries.voltage_cycles);
xlabel('Timestep');
ylabel('Voltage [V]');
title('Superimposed voltage cycles');
% %Plot the superimposed current cycles
subplot(2,1,2);
plot(output_extract_plasma_param_IV_timeseries.current_cycles);
title('Superimposed current cycles');
xlabel('Timestep');
ylabel('Current [A]');

f = figure;
%Avoid graphical errors
set(gcf,'Renderer','painters');
movegui(f,'northeast');



subplot(3,1,1)
%plot((out(1).t + timing.recMon.start)*1e3, out(1).data);
plot((out(1).t)*1e6 , out(1).data);
title('Electron Gun 1 Cathode Current');
xlabel(time_label_string);
ylabel('I_{Cathode} [A]');
%time_limits = xlim;
xlim(time_limits);

SP= 360; %your point goes here 
line([SP SP],get(gca,'YLim'),'Color',[1 0 0]);
SP= 480; %your point goes here 
line([SP SP],get(gca,'YLim'),'Color',[1 0 0]);
SP= 600; %your point goes here 
line([SP SP],get(gca,'YLim'),'Color',[1 0 0]);

% subplot(4,1,2)
% plot(time_axis,V_float,'*-k');
% title('Floating potential');
% xlabel(time_label_string);
% ylabel('V_{float} [V]');
% xlim(time_limits);

% subplot(4,1,2)
% plot(time_axis,I_sat_i_1st,'*-b');
% title('Ion saturation current (1st estimate)');
% xlabel(time_label_string);
% ylabel('I_{sat,i,1st} [A]');
% xlim(time_limits);

subplot(3,1,2)
hold on;
%plot(time_axis,T_e_max,'*-r');
plot(time_axis,T_e_min,'*-r');
%plot(time_axis,T_e_median,'*-g');
title('Electron temperature');
xlabel(time_label_string);
ylabel('T_e [eV]');
xlim(time_limits);
ylim([0 20]);

SP= 360; %your point goes here 
line([SP SP],get(gca,'YLim'),'Color',[1 0 0]);
SP= 480; %your point goes here 
line([SP SP],get(gca,'YLim'),'Color',[1 0 0]);
SP= 600; %your point goes here 
line([SP SP],get(gca,'YLim'),'Color',[1 0 0]);

subplot(3,1,3)
hold on;
%plot(time_axis,n_e_max,'*-k');
plot(time_axis,n_e_min,'*-b');
title('Electron density');
xlabel(time_label_string);
ylabel('n_e [m^{-3}]');
xlim(time_limits);

SP= 360; %your point goes here 
line([SP SP],get(gca,'YLim'),'Color',[1 0 0]);
SP= 480; %your point goes here 
line([SP SP],get(gca,'YLim'),'Color',[1 0 0]);
SP= 600; %your point goes here 
line([SP SP],get(gca,'YLim'),'Color',[1 0 0]);

