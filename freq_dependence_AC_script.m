%Similar to "extract_plasma_param_IV_timeseries", but here we do an average
%over some measurements series before fitting the data

% %TIZIANO STUFF
meas_directory = '\\pc-e5-ws-2\D\tzf\Frequency_dependence_measurements\AC_10kHz\'
voltage_sign = +1;
current_sign = -1; %(IV-Curve convention: invert current axis);
%Probe tip area
L_tip = 4.00*1e-3; %4.00 mm length
r_tip = 0.25*1e-3; %0.25 mm radius
S_tip = L_tip*pi*r_tip;
% Ion mass number
A = 39.948;%Argon

probe_tip_area = S_tip;
ion_mass_number = A;

%Select the meaurement files
m_begin = 0;
m_end   = 9;

L = m_end - m_begin + 1;

%Produce the file_names array
full_file_path_array = cell(L,1);
for n = 1 : L
    m_current = n + m_begin - 1;
    full_file_path_array{n} = strcat(meas_directory,'meas_',num2str(m_current,'%04.f'),'.h5');
end

N_samples = 5003;

data_cell_array = cell(L,1);
voltage_array = nan(L,N_samples);
current_array = nan(L,N_samples);

%Read data
for n = 1:L
    data_cell_array{n} = h5read(full_file_path_array{n},'/dataset');
    data = data_cell_array{n};
    data = data';
    voltage_array(n,:) = voltage_sign.*data(:,1);
    current_array(n,:) = current_sign.*data(:,2);
    if n == 1
        %Check metadata
        metadata = h5info(full_file_path_array{n});
        %Extract timestep from metadata [s]
        time_step = metadata.Groups.Attributes(2,1).Value(1);
        %Extract date and time from metadata
        date = metadata.Datasets.Attributes(1,1).Value;
        time = metadata.Datasets.Attributes(2,1).Value;
        %Extract uncertainty coming from ADC discretization;
        delta_voltage = signal_to_step(data(:,1));
        delta_current = signal_to_step(data(:,2));
        %Calculate sampling frequency [Hz]
        sampling_freq = 1/time_step;
    end
end

voltage = mean(voltage_array);
current = mean(current_array);


%DAQ delay in s
DAQ_delay = 20e-3;



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
        fprintf('extract_plasma_parameter_IV_timeseries: Plasma parameter extraction failed\n');
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

%Plot everything
figure;
set(gcf,'Renderer','painters');

for k = 1:NC_tot
    plot(V_axis{k},I_of_V_poly{k},'.-');
    xlim([-80 +80]);
    ylim([-0.1 3]);
    drawnow;
    pause(0.2);
end