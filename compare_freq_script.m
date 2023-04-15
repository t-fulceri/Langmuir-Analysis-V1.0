main_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\freq_dependence\';

load([main_folder, '\DAQ_delay']);
load([main_folder, '\time_step']);
load([main_folder, '\total_time_axis']);

load([main_folder, '\DC\voltage_axis']);
load([main_folder, '\DC\current_2D_array']);

load([main_folder, '\AC_10kHz\V_axis']);
load([main_folder, '\AC_10kHz\I_of_V']);
load([main_folder, '\AC_10kHz\I_of_V_poly']);
load([main_folder, '\AC_10kHz\ts_per_cycle']);
load([main_folder, '\AC_10kHz\time_offset']);
V_axis_10k = V_axis;
I_of_V_10k = I_of_V;
I_of_V_poly_10k = I_of_V_poly;
ts_per_cycle_10k = ts_per_cycle;
time_offset_10k = time_offset;

load([main_folder, '\AC_25kHz\V_axis']);
load([main_folder, '\AC_25kHz\I_of_V']);
load([main_folder, '\AC_25kHz\I_of_V_poly']);
load([main_folder, '\AC_25kHz\ts_per_cycle']);
load([main_folder, '\AC_25kHz\time_offset']);
V_axis_25k = V_axis;
I_of_V_25k = I_of_V;
I_of_V_poly_25k = I_of_V_poly;
close all;

ts_per_cycle_25k = ts_per_cycle;
time_offset_25k = time_offset;

load([main_folder, '\AC_50kHz\V_axis']);
load([main_folder, '\AC_50kHz\I_of_V']);
load([main_folder, '\AC_50kHz\I_of_V_poly']);
load([main_folder, '\AC_50kHz\ts_per_cycle']);
load([main_folder, '\AC_50kHz\time_offset']);
V_axis_50k = V_axis;
I_of_V_50k = I_of_V;
I_of_V_poly_50k = I_of_V_poly;
ts_per_cycle_50k = ts_per_cycle;
time_offset_50k = time_offset;

clear V_axis I_of_V I_of_V_poly ts_per_cycle time_offset

%PLOT EVERYTHING
f = figure;
movegui(f,'center');
set(gcf,'Renderer','painters');

%Record video
ts_per_frame = 10;
L = size(current_2D_array,3);
frames = cell(floor(L/ts_per_frame),1);
conversion_factor = 1e6; %seconds to microseconds

for ts = 1:ts_per_frame:L;
    
    subplot(2,1,1);
    time_after_daq_trigger = (ts*time_step)*conversion_factor;
    
    %DC Measurement
    plot(voltage_axis,current_2D_array(:,ts),'.-k');
    title(strcat('IV-characteristic @ time = ',num2str(time_after_daq_trigger),' µs'));
    hold on
    legend('DC measurement');
    
    %AC 10 kHz
    NC_tot_10k = length(V_axis_10k);
    for s = 1:NC_tot_10k
        if ts > time_offset_10k + (s-1)*ts_per_cycle_10k && ts <= time_offset_10k + s*ts_per_cycle_10k
            plot(V_axis_10k{s},I_of_V_10k{s},'.b');
            plot(V_axis_10k{s},I_of_V_poly_10k{s},'-b');
        end
    end
    
    %AC 25 kHz
    NC_tot_25k = length(V_axis_25k);
    for s = 1:NC_tot_25k
        if ts > time_offset_25k + (s-1)*ts_per_cycle_25k && ts <= time_offset_25k + s*ts_per_cycle_25k
            plot(V_axis_25k{s},I_of_V_25k{s},'.g');
            plot(V_axis_25k{s},I_of_V_poly_25k{s},'-g');
        end
    end
    
    %AC 50 kHz
    NC_tot_50k = length(V_axis_50k);
    for s = 1:NC_tot_50k
        if ts > time_offset_50k + (s-1)*ts_per_cycle_50k && ts <= time_offset_50k + s*ts_per_cycle_50k
            plot(V_axis_50k{s},I_of_V_50k{s},'.r');
            plot(V_axis_50k{s},I_of_V_poly_50k{s},'-r');
        end
    end
    
    
    xlabel('Voltage [V]');
    ylabel('Current [A]');
    xlim([-80 +80]);
    ylim([-0.1 3]);
    
    hold off
    
    subplot(2,1,2);
    xlim([1 L]);
    ylim([0 1]);
    
    %10 kHz snapshots
    
    drawnow;
    frames{ceil(ts/ts_per_frame)} = getframe(f);
    
end

%Save video to file
frames = cell2mat(frames);
video_directory = 'C:\Users\tzf\Desktop\VINETA II Lab\freq_dependence\';
v = VideoWriter([video_directory,'freq_dependence_video.avi'],'Uncompressed AVI');

open(v);
writeVideo(v,frames);
close(v);
