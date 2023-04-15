clear all
close all
%Access grid scan directory

% meas_directory = '\\\pc-e5-ws-2\D\tzf\Plasma_Profile_010kHz_4V_1rep\'
% %Set the number of repetition per position
% rep_per_point = 1;

% meas_directory = '\\pc-e5-ws-2\D\tzf\Plasma_Profile_050kHz_4V_5_rep_CentralZone_REPRISE\'
% meas_directory = '\\pc-e5-ws-2\D\tzf\Plasma_Profile_050kHz_4V_5_rep_CentralZone_X-Drive_ON\'
meas_directory = '\\pc-e5-ws-2\D\tzf\2018_04_19_Plasma_Profile_050kHz_4V_5_rep_CompleteProfile_X-Drive_OFF\'
%Set the number of repetition per position
rep_per_point = 5;

map_file = 'map.map'
voltage_sign = +1;
current_sign = -1;
full_file_path = strcat(meas_directory,map_file);

%Raster step
raster_step = 10;

%Read the map file
B = tdfread(full_file_path);
ind = B.x0x25_idx;
x_pos = B.x;
y_pos = B.y;
ind = ind(:,1);
x_pos = x_pos(:,1);
y_pos = y_pos(:,1);
%Number of points
NX = (max(x_pos) - min(x_pos))/raster_step + 1;
NY = (max(y_pos) - min(y_pos))/raster_step + 1;

%Organize data in multidimensional array
%Measurement number 2Dspace*1Drepetition array
meas_number_STR = NaN(NX,NY,rep_per_point);

for k = 0:length(ind)-1
    k_x = floor(k/NY) + 1;
    k_y = mod(k,NY)   + 1;

    for repetition = 1:rep_per_point
        %file_name = strcat('meas_',num2str(k+(repetition-1),'%04.f'),'.h5');
        %file_names_STR(k_x,k_y,repetition) = file_name;
        meas_number_STR(k_x,k_y,repetition) = rep_per_point*k+(repetition-1);
    end

end

%Select sub-grid to process
NX_begin = 1;
NX_end = NX;
nX = NX_end - NX_begin + 1;

NY_begin = 1;
NY_end = NY;
nY = NY_end - NY_begin + 1;

rep_begin = 1;
rep_end = 5;
nr = rep_end - rep_begin +1;

%Electron temperature 2Dspace*1Dtime*1Drepetition array
T_e_XYTR = NaN(nX,nY,100,nr);
%Electron density 2Dspace*1Dtime*1Drepetition array
n_e_XYTR = NaN(nX,nY,100,nr);
%Plasma potential 2Dspace*1Dtime*1Drepetition array
V_plasma_XYTR = NaN(nX,nY,100,nr);

%Probe tip area
L_tip = 4.00*1e-3; %4.00 mm length
r_tip = 0.25*1e-3; %0.25 mm radius
%Result:
S_tip = L_tip*pi*r_tip;

% %Ion mass number
A = 39.948;%Argon

for k_x = 1:nX
    for k_y = 1:nY
        for k_r = 1:nr
            file_number = meas_number_STR(k_x+NX_begin-1,k_y+NY_begin-1,k_r + rep_begin-1);
            filename = strcat('meas_',num2str(file_number,'%04.f'),'.h5')
            [ return_status, output_extract_plasma_param_IV_timeseries ] = extract_plasma_param_IV_timeseries( meas_directory, filename, voltage_sign, current_sign, S_tip, A);
            NC_tot = output_extract_plasma_param_IV_timeseries.NC_tot;
            T_e_XYTR(k_x,k_y,1:NC_tot,k_r) = output_extract_plasma_param_IV_timeseries.T_e_max;
            n_e_XYTR(k_x,k_y,1:NC_tot,k_r) = output_extract_plasma_param_IV_timeseries.n_e_max;
            
            T_e_XYTR = T_e_XYTR(:,:,1:NC_tot,:);
            n_e_XYTR = n_e_XYTR(:,:,1:NC_tot,:);
        end
    end
end

T_e_XYT = nanmean(T_e_XYTR,4);
n_e_XYT = nanmean(n_e_XYTR,4);