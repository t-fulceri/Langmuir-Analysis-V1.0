close all;

%DC MEASUREMENTS-----------------------------------------------------------

%Select directory
meas_directory = '\\pc-e5-ws-2\D\tzf\Frequency_dependence_measurements\DC_second_run\';
voltage_sign = +1;
current_sign = -1; %(IV-Curve convention: invert current axis);

%Select the meaurement files
m_begin = 0;
m_end   = 149;

L = m_end - m_begin + 1;

%Produce the file_names array
full_file_path_array = cell(L,1);
for ts_per_frame = 1 : L
    m_current = ts_per_frame + m_begin - 1;
    full_file_path_array{ts_per_frame} = strcat(meas_directory,'meas_',num2str(m_current,'%04.f'),'.h5');
end

%Total number of voltage levels
N_vl = 15;
%Number of measurements per voltage level
M_pvl = 10;
%Produce the file_names matrix
full_file_path_matrix = (reshape(full_file_path_array,M_pvl,N_vl))';

clear full_file_path_array;

length_taken = false;
%Fill the matrices with the proper values
for i = 1:N_vl
    for j = 1:M_pvl
        %Read data
        data = h5read(full_file_path_matrix{i,j},'/dataset');
        if ~length_taken
            L = length(data);
            voltage_3D_array = nan(N_vl,M_pvl,L);
            current_3D_array = nan(N_vl,M_pvl,L);
            length_taken = true;
        end
        %Transpose
        data = data';
        %Extract voltage fixed mean value
        voltage_3D_array(i,j,:) = data(:,1);
        current_3D_array(i,j,:) = data(:,2);
    end
end

current_3D_array = current_sign*current_3D_array;

voltage_2D_array = mean(voltage_3D_array,2);
current_2D_array = mean(current_3D_array,2);

voltage_axis = mean(voltage_2D_array,3);

%PLOT EVERYTHING
f = figure;
movegui(f,'center');
set(gcf,'Renderer','painters');

%Record video
ts_per_frame = 10;
frames = cell(floor(L/ts_per_frame),1);

for ts = 1:ts_per_frame:L;
    plot(voltage_axis,current_2D_array(:,ts),'.-');
    title(strcat('IV-characteristic @ time\_step = ',num2str(ts)));
    xlabel('Voltage [V]');
    ylabel('Current [A]');
    xlim([-80 +80]);
    ylim([-0.1 3]);
    drawnow;
    frames{ceil(ts/ts_per_frame)} = getframe(f);
end

%Save video to file
frames = cell2mat(frames);
video_directory = 'C:\Users\tzf\Desktop\VINETA II Lab\freq_dependence\'
v = VideoWriter([video_directory,'DC_meas_video.avi'],'Uncompressed AVI');

open(v);
writeVideo(v,frames);
close(v);

%--------------------------------------------------------------------------


