%Clear workspace
clear all;
%Close all figures
close all;

%TIZIANO STUFF
meas_directory = '\\pc-e5-ws-2\D\tzf\Plasma_Profile_050kHz_4V_5_rep_CentralZone_REPRISE\'
filename = 'meas_0303.h5';
voltage_sign = +1;
current_sign = -1; %(IV-Curve convention: invert current axis);
L_tip = 4.00*1e-3;
r_tip = 0.25*1e-3;
%Probe tip surface area
S_tip = L_tip*2*pi*r_tip;
%Ion mass number (depends on the gas used)
A = 39.948; %ARGON

full_file_path = strcat(meas_directory,filename);

%Check metadata
metadata = h5info(full_file_path);
%Extract timestep from metadata [s]
time_step = metadata.Groups.Attributes(2,1).Value(1);
%Extract date and time from metadata
date = metadata.Datasets.Attributes(1,1).Value
time = metadata.Datasets.Attributes(2,1).Value

%Read data
data = h5read(full_file_path,'/dataset');
%Transpose
data = data';
%Extract voltage
voltage = data(:,1);
%Optional inversion of voltage signal (depends on the data source)
voltage = voltage_sign.*voltage;
%Extract current
current = data(:,2);
%Optional inversion of current signal (depends on the data source)
current = current_sign.*current;

%Remove parasitic
%current----------------------------------------------------
parasitic_meas_directory = '\\pc-e5-ws-2\D\tzf\Vacuum_LC measurements\Frequency_Sweep_010-100kHz_50mV_CurrentMonitor_Resolution\';
%Select appropriate parasitic measurement based on frequency:
%50 kHz case
TF = strfind(meas_directory,'050kHz');
if ~isempty(TF)
    parasitic_filename = 'meas_0004.h5';
else
    %100 kHz case
    TF = strfind(meas_directory,'100kHz');
    if ~isempty(TF)
        parasitic_filename = 'meas_0009.h5';
    end
end

if isempty(TF)
    %No parasitic measurement exist for the selected bias frequency
    fprintf('No parasitic measurement exist for the selected bias frequency\n');
    fitted_parasitic = zeros(1,length(current));
else
    [ fitted_parasitic ] = extract_parasitic_current( parasitic_meas_directory, parasitic_filename , current_sign);
end
%Subtract the extracted parasitic current
current = current - fitted_parasitic;
%-------------------------------------------------------------------------------

f = figure;
movegui(f,'northwest');
set(gcf,'Renderer','painters');
subplot(2,1,1);
plot(voltage);
title('Voltage signal');
subplot(2,1,2);
plot(current + fitted_parasitic);
title('Current signal');
hold on;
plot(current);
plot(fitted_parasitic);


sampling_freq = 1/time_step;
[voltage_cycles]             =  timeseries_split_sinefit(voltage,voltage,sampling_freq);
[current_cycles,time_offset] =  timeseries_split_sinefit(voltage,current,sampling_freq);

NC_tot = size(voltage_cycles,2);

k = 10

half_cycle_length = floor(length(voltage_cycles)/2);
voltage_cycles_rise = voltage_cycles(1:half_cycle_length,:);
current_cycles_rise = current_cycles(1:half_cycle_length,:);
voltage_cycles_fall = voltage_cycles(half_cycle_length+1:end,:);
current_cycles_fall = current_cycles(half_cycle_length+1:end,:);

f = figure;
movegui(f,'north');
set(gcf,'Renderer','painters');
subplot(2,1,1);
plot(voltage_cycles(:,k),'-');
title(['Voltage cycle number ' num2str(k)]);
subplot(2,1,2);
plot(current_cycles(:,k),'-');
title(['Current cycle number ' num2str(k)]);

[ output_args ] = reconstruct_curve_improved( voltage_cycles_rise(:,k), current_cycles_rise(:,k));
V_rise = output_args.x;
I_rise = output_args.y_of_x;

[ output_args ] = reconstruct_curve_improved( voltage_cycles_fall(:,k), current_cycles_fall(:,k));
V_fall = output_args.x;
I_fall = output_args.y_of_x;

f = figure;
movegui(f,'northeast');
set(gcf,'Renderer','painters');
plot(V_rise,I_rise,'.-r');
hold on;
plot(V_fall,I_fall,'.-k');