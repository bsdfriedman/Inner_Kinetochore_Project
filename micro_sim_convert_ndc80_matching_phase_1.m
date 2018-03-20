function micro_sim_convert_ndc80_matching_phase_1(txt_file_dir)

%% Identify the files
% find the .txt files you're looking for
cd(txt_file_dir);
% identify the timepoint file to get the structure name
timepoint_file = dir('*timepoints.txt');
timepoint_file = timepoint_file.name;
% split the timepoint file
file_struc = strsplit(timepoint_file,'_');
% identify the color file names
color_struc = sprintf('ndc80_circles_%d_%d_%d_colors_r*',str2double(file_struc(3)),str2double(file_struc(4)),str2double(file_struc(5)));
color_dir = dir(color_struc);
color_rad = zeros([size(color_dir,1) 1]); % preallocate
% pull out all the radii numbers
for z = 1:size(color_dir,1)
    rad_split = strsplit(color_dir(z).name,{'_r','.'});
    color_rad(z,1) = str2double(rad_split(2));
end
% sort the radii numbers
color_rad = sort(color_rad);
 
%% Run the files through the simulator
% make the images from the cen color file
for z = 1:size(color_rad,1)
    system(sprintf('python.exe ParseBrownian.py -PSF PSF_gain_GFP.txt -out XML_ndc80_circ_r%d -width 100 -height 100 -pixel_size 64 -voxel_depth 300 -focal_planes 7 ndc80_circles_%d_%d_%d_colors_r%d.txt ndc80_circles_%d_%d_%d_timepoints.txt',color_rad(z,1),str2double(file_struc(3)),str2double(file_struc(4)),str2double(file_struc(5)),color_rad(z,1),str2double(file_struc(3)),str2double(file_struc(4)),str2double(file_struc(5))));
    system(sprintf('python.exe BrownianXMLtoTIFF.py -green -out tiff_ndc80_circ_r%d XML_ndc80_circ_r%d',color_rad(z,1),color_rad(z,1)));
end

%% Run the files through the simulator in -yzx
% make the images from the cen color file
for z = 1:size(color_rad,1)
    system(sprintf('python.exe ParseBrownian.py -PSF PSF_gain_GFP.txt -out XML_ndc80_circ_r%d_yz -width 100 -height 100 -pixel_size 64 -voxel_depth 300 -focal_planes 7 -yzx ndc80_circles_%d_%d_%d_colors_r%d.txt ndc80_circles_%d_%d_%d_timepoints.txt',color_rad(z,1),str2double(file_struc(3)),str2double(file_struc(4)),str2double(file_struc(5)),color_rad(z,1),str2double(file_struc(3)),str2double(file_struc(4)),str2double(file_struc(5))));
    system(sprintf('python.exe BrownianXMLtoTIFF.py -green -out tiff_ndc80_circ_r%d_yz XML_ndc80_circ_r%d_yz',color_rad(z,1),color_rad(z,1)));
end

