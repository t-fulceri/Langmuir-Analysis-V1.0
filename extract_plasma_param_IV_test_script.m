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

k = 14

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

%"RISE" PART
[ return_status, output_extract_plasma_param_IV ] = extract_plasma_param_IV( voltage_cycles_rise(:,k), current_cycles_rise(:,k), S_tip, A)
if return_status == 0
    V_axis = output_extract_plasma_param_IV.V_axis;
    I_of_V = output_extract_plasma_param_IV.I_of_V;
    delta_I = output_extract_plasma_param_IV.delta_I;
    I_of_V_poly = output_extract_plasma_param_IV.I_of_V_poly;
    ind_restricted = output_extract_plasma_param_IV.ind_restricted;
    shift = output_extract_plasma_param_IV.shift;
    log_I_of_V_poly = output_extract_plasma_param_IV.log_I_of_V_poly;
    D_log_I_of_V_poly = output_extract_plasma_param_IV.D_log_I_of_V_poly;
    D2_log_I_of_V_poly = output_extract_plasma_param_IV.D2_log_I_of_V_poly;
    
    T_e_lin_fit = output_extract_plasma_param_IV.T_e_lin_fit
    C_lin_fit   = output_extract_plasma_param_IV.C_lin_fit
    ind_exp_cell= output_extract_plasma_param_IV.ind_exp_cell;
    fitresult   = output_extract_plasma_param_IV.fitresult;
    gof         = output_extract_plasma_param_IV.fitresult;

    %Normalize it
    N_V = max(V_axis) - min(V_axis);
    N_I = max(I_of_V) - min(I_of_V);
    
    if length(T_e_lin_fit) > 1
%         fprintf('Warning: undecidable electron temperature fit!\n');
%         fprintf('Second value has been selected\n');
%         i = 2
%         T_e_lin_fit = T_e_lin_fit(i);
%         C_lin_fit = C_lin_fit(i);
%         ind_exp_cell = ind_exp_cell(i);

        fprintf('More than one candidate fit.\n');
        fprintf('Selecting the best one based goodness parameter\n');
        goodness = nan(length(T_e_lin_fit),1);
        for i = 1:length(T_e_lin_fit)
            goodness(i) = output_extract_plasma_param_IV.gof{i}.rmse;
        end
        [m,m_l] = min(goodness);
        
        T_e_lin_fit = T_e_lin_fit(m_l);
        C_lin_fit = C_lin_fit(m_l);
        ind_exp_cell = ind_exp_cell(m_l);
    end
    
    ind_exp = ind_exp_cell{1};
    
    V_tf = V_axis(ind_exp);
    %I_tf = I_of_V(ind_exp);
    I_tf = I_of_V_poly(ind_exp);
    
    V_tf_norm = V_tf/N_V;
    I_tf_norm = I_tf/N_I;

    %PLOT EVERYTHING
    f = figure;
    movegui(f,'northeast');
    set(gcf,'Renderer','painters');
    errorbar(V_axis,I_of_V,delta_I,'r.');
    hold on;
    plot(V_axis(ind_exp),I_of_V(ind_exp),'or');
    plot(V_axis,I_of_V_poly,'-r');
    voltage_limits = xlim;
    current_limits = ylim;
    xlabel('Voltage');
    ylabel('Current');
    title(['I(V) curve, and fit of its first part, (cycle number ' num2str(k) ')']);
    
    hold on;
    
%     f = figure;
%     movegui(f,'south');
%     set(gcf,'Renderer','painters');
%     
%     subplot(3,1,1);
%     plot(V_axis(ind_restricted),log_I_of_V_poly(ind_restricted),'-');
%     hold on;
%     fitted_log_I_of_V = C_lin_fit + (1/T_e_lin_fit)*V_axis;
%     plot_range = ind_exp(1)-30:ind_exp(end)+30;
%     plot(V_axis(ind_exp),log_I_of_V_poly(ind_exp),'ok');
%     plot(V_axis(plot_range),fitted_log_I_of_V(plot_range),'-');  
%     xlabel('Voltage');
%     ylabel('ln(I)');
%     title(['log(I(V)) (approximated) (cycle number ' num2str(k) ')']);
%     xlim(voltage_limits);
%     ylim([-5 0]);
%     
%     subplot(3,1,2);
%     plot(V_axis(ind_restricted),D_log_I_of_V_poly(ind_restricted),'-');
%     hold on;
%     plot(V_axis(ind_exp),D_log_I_of_V_poly(ind_exp),'ok');
%     xlabel('Voltage');
%     ylabel('d(ln(I))/dV');
%     title(['1st Derivative of log(I(V)) (approximated) (cycle number ' num2str(k) ')']);
%     xlim(voltage_limits);
%     
%     subplot(3,1,3);
%     plot(V_axis(ind_restricted),D2_log_I_of_V_poly(ind_restricted),'-');
%     hold on;
%     plot(V_axis(ind_exp),D2_log_I_of_V_poly(ind_exp),'ok');
%     xlabel('Voltage');
%     ylabel('d^2(ln(I))/dV^2');
%     title(['2nd Derivative of log(I(V)) (approximated) (cycle number ' num2str(k) ')']);
%     xlim(voltage_limits);
    

%     f = figure;
%     movegui(f,'southeast');
%     set(gcf,'Renderer','painters');
%     plot(V_axis(ind_restricted),D_I_of_V_poly(ind_restricted),'-');
%     xlabel('Voltage');
%     ylabel('Current/Voltage');
%     title(['dI/dV curve (approximated) (cycle number ' num2str(k) ')']);
%     xlim(time_limits);
    
%     f = figure;
%     movegui(f,'southeast');
%     set(gcf,'Renderer','painters');
%     errorbar(V_axis(exp_ind),I_of_V(exp_ind),delta_I(exp_ind),'.');
%     hold on;
%     plot(V_axis(exp_ind),fitted_I_of_V2(exp_ind),'-');
%     xlabel('Voltage');
%     ylabel('Current');
%     xlim(voltage_limits);
%     ylim(current_limits);



end

%"FALL" PART
[ return_status, output_extract_plasma_param_IV ] = extract_plasma_param_IV( voltage_cycles_fall(:,k), current_cycles_fall(:,k), S_tip, A)
if return_status == 0
    V_axis = output_extract_plasma_param_IV.V_axis;
    I_of_V = output_extract_plasma_param_IV.I_of_V;
    delta_I = output_extract_plasma_param_IV.delta_I;
    I_of_V_poly = output_extract_plasma_param_IV.I_of_V_poly;
    ind_restricted = output_extract_plasma_param_IV.ind_restricted;
    shift = output_extract_plasma_param_IV.shift;
    log_I_of_V_poly = output_extract_plasma_param_IV.log_I_of_V_poly;
    D_log_I_of_V_poly = output_extract_plasma_param_IV.D_log_I_of_V_poly;
    D2_log_I_of_V_poly = output_extract_plasma_param_IV.D2_log_I_of_V_poly;
    
    T_e_lin_fit = output_extract_plasma_param_IV.T_e_lin_fit
    C_lin_fit   = output_extract_plasma_param_IV.C_lin_fit
    ind_exp_cell= output_extract_plasma_param_IV.ind_exp_cell;
    fitresult   = output_extract_plasma_param_IV.fitresult;
    gof         = output_extract_plasma_param_IV.fitresult;

    %Normalize it
    N_V = max(V_axis) - min(V_axis);
    N_I = max(I_of_V) - min(I_of_V);
    
    if length(T_e_lin_fit) > 1
%         fprintf('Warning: undecidable electron temperature fit!\n');
%         fprintf('Second value has been selected\n');
%         i = 2
%         T_e_lin_fit = T_e_lin_fit(i);
%         C_lin_fit = C_lin_fit(i);
%         ind_exp_cell = ind_exp_cell(i);

        fprintf('More than one candidate fit.\n');
        fprintf('Selecting the best one based goodness parameter\n');
        goodness = nan(length(T_e_lin_fit),1);
        for i = 1:length(T_e_lin_fit)
            goodness(i) = output_extract_plasma_param_IV.gof{i}.rmse;
        end
        [m,m_l] = min(goodness);
        
        T_e_lin_fit = T_e_lin_fit(m_l);
        C_lin_fit = C_lin_fit(m_l);
        ind_exp_cell = ind_exp_cell(m_l);
    end
    
    ind_exp = ind_exp_cell{1};
    
    V_tf = V_axis(ind_exp);
    %I_tf = I_of_V(ind_exp);
    I_tf = I_of_V_poly(ind_exp);
    
    V_tf_norm = V_tf/N_V;
    I_tf_norm = I_tf/N_I;

    %PLOT EVERYTHING
%     f = figure;
%     movegui(f,'southeast');
%     set(gcf,'Renderer','painters');
    errorbar(V_axis,I_of_V,delta_I,'.k');
    hold on;
    plot(V_axis(ind_exp),I_of_V(ind_exp),'ok');
    plot(V_axis,I_of_V_poly,'-k');
    voltage_limits = xlim;
    current_limits = ylim;
    xlabel('Voltage');
    ylabel('Current');
    title(['I(V) curve, and fit of its first part, (cycle number ' num2str(k) ')']);
    
%     f = figure;
%     movegui(f,'south');
%     set(gcf,'Renderer','painters');
%     
%     subplot(3,1,1);
%     plot(V_axis(ind_restricted),log_I_of_V_poly(ind_restricted),'-');
%     hold on;
%     fitted_log_I_of_V = C_lin_fit + (1/T_e_lin_fit)*V_axis;
%     plot_range = ind_exp(1)-30:ind_exp(end)+30;
%     plot(V_axis(ind_exp),log_I_of_V_poly(ind_exp),'ok');
%     plot(V_axis(plot_range),fitted_log_I_of_V(plot_range),'-');  
%     xlabel('Voltage');
%     ylabel('ln(I)');
%     title(['log(I(V)) (approximated) (cycle number ' num2str(k) ')']);
%     xlim(voltage_limits);
%     ylim([-5 0]);
%     
%     subplot(3,1,2);
%     plot(V_axis(ind_restricted),D_log_I_of_V_poly(ind_restricted),'-');
%     hold on;
%     plot(V_axis(ind_exp),D_log_I_of_V_poly(ind_exp),'ok');
%     xlabel('Voltage');
%     ylabel('d(ln(I))/dV');
%     title(['1st Derivative of log(I(V)) (approximated) (cycle number ' num2str(k) ')']);
%     xlim(voltage_limits);
%     
%     subplot(3,1,3);
%     plot(V_axis(ind_restricted),D2_log_I_of_V_poly(ind_restricted),'-');
%     hold on;
%     plot(V_axis(ind_exp),D2_log_I_of_V_poly(ind_exp),'ok');
%     xlabel('Voltage');
%     ylabel('d^2(ln(I))/dV^2');
%     title(['2nd Derivative of log(I(V)) (approximated) (cycle number ' num2str(k) ')']);
%     xlim(voltage_limits);
    

%     f = figure;
%     movegui(f,'southeast');
%     set(gcf,'Renderer','painters');
%     plot(V_axis(ind_restricted),D_I_of_V_poly(ind_restricted),'-');
%     xlabel('Voltage');
%     ylabel('Current/Voltage');
%     title(['dI/dV curve (approximated) (cycle number ' num2str(k) ')']);
%     xlim(time_limits);
    
%     f = figure;
%     movegui(f,'southeast');
%     set(gcf,'Renderer','painters');
%     errorbar(V_axis(exp_ind),I_of_V(exp_ind),delta_I(exp_ind),'.');
%     hold on;
%     plot(V_axis(exp_ind),fitted_I_of_V2(exp_ind),'-');
%     xlabel('Voltage');
%     ylabel('Current');
%     xlim(voltage_limits);
%     ylim(current_limits);



end

