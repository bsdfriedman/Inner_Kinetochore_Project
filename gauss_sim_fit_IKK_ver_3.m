function [coords_gauss,FWHM] = gauss_sim_fit_IKK_ver_3(image, gof_thresh)
% Written by Brandon Friedman
% This code finds the new gaussian fitted coordinates and FWHM for a simulated input image
% using a linescan

%% Determining the coordinates and FWHM in X
% find mean in first dimension for the x position
image_a = mat2gray(mean(image,1))';

% fits a gaussian to the coordinates
[b,gof_b] = fit((1:size(image,2))',image_a,'gauss1');

% extracts the coordinates and logs the center of the curve
coeff_x = coeffvalues(b);

% if the curve fits well, log the new center
if gof_b.rsquare > gof_thresh
    x_gauss = coeff_x(2);
    FWHM_x = (coeff_x(3)/sqrt(2)) * sqrt(8*log(2));
else
    % if it doesn't, log a dummy variable that will be excluded
    x_gauss = -1;
    FWHM_x = -1;
end

%% Determining the coordinates and FHWM in Y

% find mean in second dimension for the y position
image_c = mat2gray(mean(image,2));

% fits a gaussian to the coordinates
[d,gof_d] = fit((1:size(image,1))',image_c,'gauss1');

% extracts the coordinates and logs the center of the curve
coeff_y = coeffvalues(d);
% if the curve fits well, log the new center
if gof_d.rsquare > gof_thresh
    y_gauss = coeff_y(2);
    FWHM_y = (coeff_y(3)/sqrt(2)) * sqrt(8*log(2));
else
    % if it doesn't, log a dummy variable that will be excluded
    y_gauss = -1;
    FWHM_y = -1;
end

%% Logging the Coordinates
% puts the data into a data structure
coords_gauss = [x_gauss, y_gauss];
FWHM = [FWHM_x, FWHM_y];