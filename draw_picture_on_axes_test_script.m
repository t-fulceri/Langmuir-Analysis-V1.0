close all;

picture_dir = 'C:\Users\tzf\Desktop\VINETA II Lab\PPT_March_2018\';
file_name = 'E_gun_CAD.bmp';
complete_picture_path = [picture_dir,file_name];
E_gun_CAD_image = imread(complete_picture_path);
alpha_data = ~E_gun_CAD_image;




f = figure;
%Avoid graphical errors
set(gcf,'Renderer','painters');
movegui(gcf,'northwest');

h = imshow(E_gun_CAD_image);
set(h, 'AlphaData', alpha_data);
axis on