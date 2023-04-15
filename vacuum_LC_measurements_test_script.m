close all;


parasitic_meas_directory = '\\pc-e5-ws-2\D\tzf\Vacuum_LC measurements\Frequency_Sweep_010-100kHz_50mV_CurrentMonitor_Resolution\';
%At 50 kHz
parasitic_filename = 'meas_0004.h5';
voltage_sign = +1;
current_sign = -1; %(IV-Curve convention: invert current axis);

full_file_path = strcat(parasitic_meas_directory,parasitic_filename);

%Check metadata
metadata = h5info(full_file_path);
%Extract timestep from metadata [s]
time_step = metadata.Groups.Attributes(2,1).Value(1);
%Extract date and time from metadata
date = metadata.Datasets.Attributes(1,1).Value;
time = metadata.Datasets.Attributes(2,1).Value;

%Read data
data = h5read(full_file_path,'/dataset');
%DAQ delay in s
DAQ_delay = 20e-3;
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

sampling_freq = 1/time_step;
[voltage_cycles]             =  timeseries_split_sinefit(voltage,voltage,sampling_freq);
[current_cycles,time_offset] =  timeseries_split_sinefit(voltage,current,sampling_freq);

%Time steps per cycle
ts_per_cycle = size(voltage_cycles,1);
%Total number of cycles
NC_tot = size(voltage_cycles,2);

%Produce the time axis
conversion_factor = 1e+6; %10e6 => seconds to microseconds
time_axis = time_step*ts_per_cycle*conversion_factor*(linspace(1,NC_tot,NC_tot) - 1/2);
%time_axis = time_axis + DAQ_delay*conversion_factor;


total_time_axis = (1/sampling_freq)*conversion_factor*linspace(0,length(voltage)-1,length(voltage));
%total_time_axis = total_time_axis + DAQ_delay*conversion_factor;

k = 5

half_cycle_length = floor(length(voltage_cycles)/2);
voltage_cycles_rise = voltage_cycles(1:half_cycle_length,:);
current_cycles_rise = current_cycles(1:half_cycle_length,:);
voltage_cycles_fall = voltage_cycles(half_cycle_length+1:end,:);
current_cycles_fall = current_cycles(half_cycle_length+1:end,:);

[ output_args ] = reconstruct_curve_improved( voltage_cycles_rise(:,k), current_cycles_rise(:,k));
V_rise = output_args.x;
I_rise = output_args.y_of_x;

[ output_args ] = reconstruct_curve_improved( voltage_cycles_fall(:,k), current_cycles_fall(:,k));
V_fall = output_args.x;
I_fall = output_args.y_of_x;


f = figure;
movegui(f,'northwest');
set(gcf,'Renderer','painters');
subplot(2,1,1);
plot(total_time_axis,voltage);
title('Voltage bias signal (80 V AC @  50 kHz)');
xlabel('Time [µs]');
ylabel('Voltage [V]');
subplot(2,1,2);
plot(total_time_axis,current*1e3);
% hold on;
% plot(total_time_axis,fitted_parasitic*1e3);
title('Probe parasitic current signal');
xlabel('Time [µs]');
ylabel('Current [mA]');


f = figure;
movegui(f,'north');
set(gcf,'Renderer','painters');
subplot(2,1,1);
plot(voltage_cycles(:,k),'-');
title(['Voltage cycle number ' num2str(k)]);
subplot(2,1,2);
plot(current_cycles(:,k),'-');
title(['Current cycle number ' num2str(k)]);


f = figure;
movegui(f,'northeast');
set(gcf,'Renderer','painters');
plot(V_rise,I_rise*1e3,'.-r');
hold on;
plot(V_fall,I_fall*1e3,'.-k');
title(['Hysteresis loop at cycle number ' num2str(k)]);
xlabel('Voltage [V]');
ylabel('Current [mA]');