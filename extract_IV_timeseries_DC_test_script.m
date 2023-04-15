close all;

%Select directory
meas_directory = '\\pc-e5-ws-2\D\tzf\DC_measurements_central_point\';
voltage_sign = +1;
current_sign = -1; %(IV-Curve convention: invert current axis);

%Select the meaurement files
%According to personal logbook entry of 2018.3.15:
%from meas_0015 to meas_0029
m_begin = 15;
m_end   = 29;

L = m_end - m_begin + 1;

%Produce file_names
full_file_path_array = cell(L,1);
for n = 1 : L
    m_current = n + m_begin - 1;
    full_file_path_array{n} = strcat(meas_directory,'meas_',num2str(m_current,'%04.f'),'.h5');
end

%Cycle through the files and extract a constant voltage for each one of
%them
voltage_axis = nan(L,1);
current_matrix = nan(L,10000);

for n = 1 : L
%     %Check metadata
%     metadata = h5info(full_file_path_array{n});
%     %Extract timestep from metadata [s]
%     time_step = metadata.Groups.Attributes(2,1).Value(1);
%     %Extract date and time from metadata
%     date = metadata.Datasets.Attributes(1,1).Value;
%     time = metadata.Datasets.Attributes(2,1).Value;

    %Read data
    data = h5read(full_file_path_array{n},'/dataset');
    %Take note of number of samples
    n_samples = size(data,2);
    %Transpose
    data = data';
    %Extract voltage fixed mean value
    voltage_axis(n) = mean(data(:,1));
    %Optional inversion of voltage signal (depends on the data source)
    voltage_axis(n) = voltage_sign.*voltage_axis(n);
    %Extract current timeseries
    current_matrix(n,1:n_samples) = data(:,2);

end

%Optional inversion of current signal (depends on the data source)
current_matrix = current_sign.*current_matrix;
%Re-size the data matrix for current
current_matrix = current_matrix(:,1:n_samples);


%PLOT EVERYTHING
f = figure;
movegui(f,'center');
set(gcf,'Renderer','painters');

for ts = 1:n_samples;
    plot(voltage_axis,current_matrix(:,ts),'.-');
    title(strcat('IV-characteristic @ time\_step = ',num2str(ts)));
    xlabel('Voltage [V]');
    ylabel('Current [A]');
    drawnow;
end


