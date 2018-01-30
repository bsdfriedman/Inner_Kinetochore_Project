function [mask_f] = mask_maker_footprint_IKK_ver3(mask, coords, image_h, image_w)
% Written by Brandon Friedman

% invert the mask
mask_a = ~mask;

% add padding
mask_b = padarray(mask_a, [0 coords(1)-1-((size(mask,2)-1)/2)], 1, 'pre');
mask_c = padarray(mask_b, [0 image_w-coords(1)-((size(mask,2)-1)/2)], 1, 'post');
mask_d = padarray(mask_c, coords(2)-1-((size(mask,1)-1)/2), 1, 'pre');
mask_f = padarray(mask_d, image_h-coords(2)-((size(mask,1)-1)/2), 1, 'post');