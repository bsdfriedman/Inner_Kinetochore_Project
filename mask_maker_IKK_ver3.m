function [mask] = mask_maker_IKK_ver3(x_width, y_width, pad, theta)
% function to create a padded binary matrix, core is ones, padding is zeros.
core = ones([y_width,x_width]);
bin = padarray(core,[pad pad]);
 
% rotate the binary mask by theta
mask = imrotate(bin,theta,'crop');