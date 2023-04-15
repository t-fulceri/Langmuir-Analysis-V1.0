%Close all figures
close all;


%Add directories for functions
addpath('\\pcvineta3\Data\MATLAB\files');
addpath('\\pcvineta3\Data\MATLAB\files\HDF5-tools');

%Reconnection monitor output files directory
recmon_output_dir = '\\pc-e5-ws-2\D\reconnection_monitor\';

filename = 'rcc_214137.h5';

full_file_path = strcat(recmon_output_dir,filename);

[out,timing]=readRecMon(full_file_path, 'rogCalib.txt');

f = figure;
movegui(f,'north');
set(gcf,'Renderer','painters');
plot(out(1).t + timing.recMon.start, out(1).data);

timing.triggers.time(4)