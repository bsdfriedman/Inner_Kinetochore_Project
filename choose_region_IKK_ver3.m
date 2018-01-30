function [region] = choose_region_IKK_ver3(filename,x_range,y_range,sample_plane,num_z,padding)
% Written by Brandon Friedman
% this program crops the plane to select the appropriate cell region

% get the file info
info = imfinfo(filename);

% choose the plane regions
plane_ceil = ceil(sample_plane/num_z)*num_z;
plane = zeros([max([info.Height]) max([info.Width]) num_z]);
plane_pad = zeros([max([info.Height])+(2*padding) max([info.Width])+(2*padding) num_z]);

% set the plane counter to loop through planes
plane_counter = 1;

% open up the planes
for n = plane_ceil-num_z+1:plane_ceil
    plane(:,:,plane_counter)= imread(filename, n, 'Info', info);
    plane_counter = plane_counter + 1;
end

% add padding to the planes
for n = 1:num_z
    plane_pad(:,:,n) = padarray(plane(:,:,n),[padding padding]);
end

% select the region you want
region = plane_pad(y_range(1):y_range(2),x_range(1):x_range(2),:);