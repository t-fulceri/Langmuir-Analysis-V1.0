function [ return_status, output_args ] = extract_plasma_param_IV_timeseries( meas_directory, filename, voltage_sign, current_sign, probe_tip_area, ion_mass_number)
%ANALYZE IV TIMESERIES Extracts plasma parameters as a function of time
%from two paired periodic timeseries signals (voltage and current)

full_file_path = strcat(meas_directory,filename);

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
%Extract uncertainty coming from ADC discretization;
delta_voltage = signal_to_step(voltage);
delta_current = signal_to_step(current);

%Calculate sampling frequency [Hz]
sampling_freq = 1/time_step;


%Optional parasitic current removal
remove_parasitic = false;

if remove_parasitic
    fprintf('Parasitic current subtraction is ON\n');
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
        fprintf('Parasitic current has been subtracted\n');
        [ fitted_parasitic ] = extract_parasitic_current( parasitic_meas_directory, parasitic_filename , current_sign );
    end
    %Subtract the extracted parasitic current
    current = current - fitted_parasitic;
    %-------------------------------------------------------------------------------
else
    fprintf('Parasitic current subtraction is OFF\n');
end

%Split the timeseries signal into chunks of equal duration
%(each corresponding to one voltage cycle)
[voltage_cycles]             =  timeseries_split_sinefit(voltage,voltage,sampling_freq);
[current_cycles,time_offset] =  timeseries_split_sinefit(voltage,current,sampling_freq);

%Time steps per cycle
ts_per_cycle = size(voltage_cycles,1);
%Total number of cycles
NC_tot = size(voltage_cycles,2);

%Produce the time axis
conversion_factor = 1e+6; %10^6 => seconds to microseconds
time_axis = time_step*ts_per_cycle*conversion_factor*(linspace(1,NC_tot,NC_tot) - 1/2);
time_axis = time_axis + DAQ_delay*conversion_factor;


total_time_axis = (1/sampling_freq)*conversion_factor*linspace(0,length(voltage)-1,length(voltage));
total_time_axis = total_time_axis + DAQ_delay*conversion_factor;

%NEW METHOD (2018-01-30)
%Just take the half of a cycle length
%(which was calculated through fft method)
half_cycle_length = floor(length(voltage_cycles)/2);

voltage_cycles_rise = voltage_cycles(1:half_cycle_length,:);
current_cycles_rise = current_cycles(1:half_cycle_length,:);
voltage_cycles_fall = voltage_cycles(half_cycle_length+1:end,:);
current_cycles_fall = current_cycles(half_cycle_length+1:end,:);


%Prepare output variables for IV axes
V_axis = cell(NC_tot,1);
I_of_V = cell(NC_tot,1);
delta_I = cell(NC_tot,1);
I_of_V_poly = cell(NC_tot,1);


%FLOATING POTENTIAL
V_float = nan(NC_tot,1);
%ION SATURATION CURRENT (1st estimate)
I_sat_i_1st = nan(NC_tot,1);
%ELECTRON TEMPERATURE
T_e_lin_fit = cell(NC_tot,1);
%ELECTRON DENSITY
n_e_derived = cell(NC_tot,1);

%By default, the function is returning the "something wrong" state
return_status = 1;

for k = 1:NC_tot
    fprintf('Cycle_number = %d\n',k);
     try
        [ extract_plasma_param_IV_return_status, output_extract_plasma_param_IV ] = extract_plasma_param_IV( voltage_cycles_rise(:,k), current_cycles_rise(:,k), probe_tip_area, ion_mass_number);
        if extract_plasma_param_IV_return_status == 0
            
            V_axis{k} = output_extract_plasma_param_IV.V_axis;
            I_of_V{k} = output_extract_plasma_param_IV.I_of_V;
            delta_I{k} = output_extract_plasma_param_IV.delta_I;
            I_of_V_poly{k} = output_extract_plasma_param_IV.I_of_V_poly;
            shift = output_extract_plasma_param_IV.shift;
            
            V_float(k) = output_extract_plasma_param_IV.V_float;
            I_sat_i_1st(k) = output_extract_plasma_param_IV.I_sat_i_1st;
            T_e_lin_fit{k} = output_extract_plasma_param_IV.T_e_lin_fit;
            n_e_derived{k} = output_extract_plasma_param_IV.n_e_derived;
            
            %Extract max min and median of T_e and n_e
            L = length(T_e_lin_fit);
            T_e_min = NaN(L,1);
            T_e_max = NaN(L,1);

            for i = 1:L
                T_arr = T_e_lin_fit{i};
                switch length(T_arr)
                    case 0
                       T_e_min(i) = nan;
                       T_e_max(i) = nan;
                    case 1
                        T_e_min(i) = T_arr(1);
                        T_e_max(i) = T_arr(1);
                    otherwise         
                        T_e_min(i) = min(T_arr);
                        T_e_max(i) = max(T_arr);
                end
            end
            T_e_median = (T_e_min + T_e_max)./2;

            L = length(n_e_derived);
            n_e_min = NaN(L,1);
            n_e_max = NaN(L,1);

            for i = 1:L
                n_arr = n_e_derived{i};
                switch length(n_arr)
                    case 0
                       n_e_min(i) = nan;
                       n_e_max(i) = nan;
                    case 1
                        n_e_min(i) = n_arr(1);
                        n_e_max(i) = n_arr(1);
                    otherwise         
                        n_e_min(i) = min(n_arr);
                        n_e_max(i) = max(n_arr);
                end
            end
            n_e_median = (T_e_min + T_e_max)./2;
            
        else
            shift = NaN;
            
            V_axis{k} = [];
            I_of_V{k} = [];
            delta_I{k} = [];
            I_of_V_poly{k} = [];
            
            V_float(k) = NaN;
            I_sat_i_1st(k) = NaN;
            T_e_lin_fit{k} = [];
            n_e_derived{k} = [];
            
            T_e_min = NaN;
            T_e_max = NaN;
            T_e_median = NaN;
            n_e_min = NaN;
            n_e_max = NaN;
            n_e_median = NaN;
            
        end
    catch
        fprintf('extract_plasma_parameter_IV_timeseries: Plasma parameter extraction extraction failed\n');
        shift = NaN;
        
        V_axis{k} = [];
        I_of_V{k} = [];
        delta_I{k} = [];
        I_of_V_poly{k} = [];
        
        V_float(k) = NaN;
        I_sat_i_1st(k) = NaN;
        T_e_lin_fit{k} = [];
        n_e_derived{k} = [];
        
        T_e_min = NaN;
        T_e_max = NaN;
        T_e_median = NaN;
        n_e_min = NaN;
        n_e_max = NaN;
        n_e_median = NaN;
    end
end

output_args.date = date;
output_args.time = time;
output_args.time_step = time_step;
output_args.sampling_freq = sampling_freq;
output_args.ts_per_cycle = ts_per_cycle;
output_args.NC_tot = NC_tot;

output_args.voltage = voltage;
output_args.current = current;
output_args.delta_voltage = delta_voltage;
output_args.delta_current = delta_current;
output_args.voltage_cycles = voltage_cycles;
output_args.current_cycles = current_cycles;

output_args.voltage_cycles_rise = voltage_cycles_rise;
output_args.current_cycles_rise = current_cycles_rise;
output_args.voltage_cycles_fall = voltage_cycles_fall;
output_args.current_cycles_fall = current_cycles_fall;

output_args.time_axis = time_axis;
output_args.total_time_axis = total_time_axis;

output_args.shift = shift;

output_args.V_axis{k} = V_axis{k};
output_args.I_of_V{k} = I_of_V{k};
output_args.delta_I{k} = delta_I{k};
output_args.I_of_V_poly{k} = I_of_V_poly{k};

output_args.V_float = V_float;
output_args.I_sat_i_1st = I_sat_i_1st;
output_args.T_e_lin_fit = T_e_lin_fit;
output_args.n_e_derived = n_e_derived;

output_args.T_e_min = T_e_min;
output_args.T_e_max = T_e_max;
output_args.T_e_median = T_e_median;
output_args.n_e_min = n_e_min;
output_args.n_e_max = n_e_max;
output_args.n_e_median = n_e_median;


%If you arrived here, it means that everything went ok...
return_status = 0;

end

