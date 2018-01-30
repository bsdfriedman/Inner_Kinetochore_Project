function [cropped_image] = crop_image_IKK_ver3(image, center, dim_y, dim_x)
% Written by Brandon Friedman

% crop a region of size 2dim+1 in y and 2dim+1 in x around the image center
cropped_image = image(center(2)-dim_y:center(2)+dim_y,center(1)-dim_x:center(1)+dim_x);

end