function [FWHM,gof,image] = FHWM_sim_images(filename)

% parameters
mask_hy = 15; % height of the mask for the Y-FWHM measurement
mask_wy = 7; % width of the mask for the Y-FWHM measurement

% open up the file
info = imfinfo(filename);
plane(:,:,1)= double(imread(filename, 4, 'Info', info));

% find the image center
[~, max_idx] = max(plane(:));
[bp_X,bp_Y]=ind2sub(size(plane),max_idx);

% crop the image
image = crop_image_IKK_ver3(plane, [bp_X bp_Y], (mask_hy-1)/2, (mask_wy-1)/2);

% log the image height
hght = size(image,1);

% % adds noise
% a = randi(120,15,7);
% image = image+a;
% imshow(image,[]);

% perform the FHWM calculation
image_a = sum(image,2);
z = (1:size(image,1))';
equ = fittype('a1*exp(-((x-b1)/c1)^2)+z');
image_b = mat2gray(image_a);
[g,gof] = fit(z,image_b,equ,'Lower',[0.7,((hght+1)/2)-(0.25*hght),1,0],'Upper',[1.3,((hght+1)/2)+(0.25*hght),hght*0.8,0.3],'StartPoint',[1,(hght+1)*0.5,hght*0.3,0]);
coeff = coeffvalues(g);
FWHM = (coeff(3)/sqrt(2)) * sqrt(8*log(2));
