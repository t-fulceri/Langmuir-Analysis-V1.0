close all;
%WARNING: MUST BE EXECUTED AFTER A SUCCESSFUL CALL TO
%extract_plasma_param_IV_timeseries

%Load T_e and n_e profiles

%profile_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\Erlangen_2018\Matlab_profiles\'
% load([profile_folder, 'n_e_min_XYT']);
% load([profile_folder, 'T_e_min_XYT']);
% T_scale = [3 10];
% n_scale = [1 5]*1e19;

profile_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\Reconstructed_Profiles\2018_03_22\';
load([profile_folder, 'n_e_XYT']);
load([profile_folder, 'T_e_XYT']);
T_scale = [0 20];
n_scale = [0 6]*1e19;

%Smooth the data a little bit...
T_e_XYT = smooth3(T_e_XYT,'gaussian',3);
n_e_XYT = smooth3(n_e_XYT,'gaussian',3);

%snapshots = [6 ,12, 18];
%snapshots = [8];
snapshots = [4 ,8, 12, 16];
Ns_tot = length(snapshots);

f = figure;
%Avoid graphical errors
set(gcf,'Renderer','painters');
movegui(gcf,'northwest');

%Load the E_gun CAD image
picture_dir = 'C:\Users\tzf\Desktop\VINETA II Lab\PPT_March_2018\';
file_name = 'E_gun_CAD.bmp';
complete_picture_path = [picture_dir,file_name];
E_gun_CAD_image = imread(complete_picture_path);
E_gun_CAD_image = E_gun_CAD_image(203:203+278,203:203+278);
scale_factor = 100/size(E_gun_CAD_image,1);
E_gun_CAD_image = imresize(E_gun_CAD_image,scale_factor);
alpha_data = ~E_gun_CAD_image;

for ns = 1:Ns_tot
    
    subplot(2,Ns_tot,ns);
    hold on;

    
    [C,h] = contourf(98:10:198,65:10:165,T_e_XYT(:,:,snapshots(ns)));
    set(h,'LineColor','none')
    
%     hp = imshow(E_gun_CAD_image,'InitialMagnification','fit');
%     set(hp, 'AlphaData', alpha_data);
    
    title(strcat('Electron Temperature profile @ t = ',num2str(time_axis(snapshots(ns))),' 탎'));
    xlabel('X position [mm]');
    ylabel('Y position [mm]');
    pbaspect([1 1 1]);
    colormap(gca,hot);
    caxis(T_scale);
    h = colorbar(gca);
    ylabel(h, 'T_e [eV]');
    

    
    hold off;
    
    subplot(2,Ns_tot,ns + Ns_tot);
    [C,h] = contourf(98:10:198,65:10:165,n_e_XYT(:,:,snapshots(ns)));
    set(h,'LineColor','none')
    title(strcat('Electron Density profile @ t = ',num2str(time_axis(snapshots(ns))),' 탎'));
    xlabel('X position [mm]');
    ylabel('Y position [mm]');
    pbaspect([1 1 1]);
    colormap(gca,parula);
    caxis(n_scale)
    h = colorbar(gca);
    ylabel(h, 'n_e [m^{-3}]');
    
end

f = figure;
%Avoid graphical errors
set(gcf,'Renderer','painters');
%movegui(gcf,'center');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);


%Record video
frames = cell(length(time_axis),1);

for ts = 1:length(time_axis)

    subplot(2,2,[1,2]);
    plot((out(1).t)*1e6 , out(1).data);
    title('Electron Gun 1 Cathode Current');
    xlabel(time_label_string);
    ylabel('I_{Cathode} [A]');
    %time_limits = xlim;
    xlim(time_limits);

    SP= time_axis(ts); %your point goes here 
    line([SP SP],get(gca,'YLim'),'Color',[1 0 0]);

    subplot(2,2,3);
    [C,h] = contourf(98:10:198,65:10:165,T_e_XYT(:,:,ts));
    set(h,'LineColor','none')
    title(strcat('Electron Temperature profile @ t = ',num2str(time_axis(ts)),' 탎'));
    xlabel('X position [mm]');
    ylabel('Y position [mm]');
    pbaspect([1 1 1]);
    colormap(gca,hot);
    caxis(T_scale);
    h = colorbar(gca);
    ylabel(h, 'T_e [eV]');

    subplot(2,2,4);
    [C,h] = contourf(98:10:198,65:10:165,n_e_XYT(:,:,ts));
    set(h,'LineColor','none')
    title(strcat('Electron Density profile @ t = ',num2str(time_axis(ts)),' 탎'));
    xlabel('X position [mm]');
    ylabel('Y position [mm]');
    pbaspect([1 1 1]);
    colormap(gca,parula);
    caxis(n_scale)
    h = colorbar(gca);
    ylabel(h, 'n_e [m^{-3}]');
    
    drawnow;
    frames{ts} = getframe(f);

end

%Save video to file
frames = cell2mat(frames);
v = VideoWriter([profile_folder,'profiles_video.avi'],'Uncompressed AVI');

open(v);
writeVideo(v,frames);
close(v);