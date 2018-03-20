function [coords_gauss] = gauss_center_calc_ndc80_matching(image, gof_thresh)
% Written by Brandon Friedman
% This code finds the new gaussian fitted coordinates for an input image
% using a linescan

%% Determining the coordinates in X
% find mean in first dimension for the x position
image_a = mat2gray(mean(image,1))';

% parameters
hght = size(image,2);

% fits a gaussian to the coordinates
equ = fittype('a1*exp(-((x-b1)/c1)^2)+z');
[b,gof_b] = fit((1:size(image,2))',image_a,equ,'Lower',[0.7,((hght+1)/2)-(0.25*hght),1,0],'Upper',[1.3,((hght+1)/2)+(0.25*hght),hght*0.8,0.3],'StartPoint',[1,(hght+1)*0.5,hght*0.3,0]);

% extracts the coordinates and logs the center of the curve
coeff_x = coeffvalues(b);
% if the curve fits well, log the new center
if gof_b.rsquare > gof_thresh
    x_gauss = coeff_x(2);
else
    % if it doesn't, log a dummy variable that will be excluded
    x_gauss = -1000000;
end

%% Determining the coordinates in Y

% find mean in second dimension for the y position
image_c = mat2gray(mean(image,2));

% fits a gaussian to the coordinates
[d,gof_d] = fit((1:size(image,1))',image_c,equ,'Lower',[0.7,((hght+1)/2)-(0.25*hght),1,0],'Upper',[1.3,((hght+1)/2)+(0.25*hght),hght*0.8,0.3],'StartPoint',[1,(hght+1)*0.5,hght*0.3,0]);

% extracts the coordinates and logs the center of the curve
coeff_y = coeffvalues(d);
% if the curve fits well, log the new center
if gof_d.rsquare > gof_thresh
    y_gauss = coeff_y(2);
else
    % if it doesn't, log a dummy variable that will be excluded
    y_gauss = -1000000;
end

%% Logging the Coordinates
% puts the data into a data structure
coords_gauss = [x_gauss, y_gauss];