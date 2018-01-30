function [FWHM] = FWHM_calculate_y_IKK_ver3(image,hght,gof_thresh)
% Written by Brandon Friedman

% find mean in second dimension
image_a = sum(image,2);

% make x axis
z = (1:size(image,1))';

% generate a gaussian curve template with a +z correction to account for noise
equ = fittype('a1*exp(-((x-b1)/c1)^2)+z');

% fit the curve to the gaussian
image_b = mat2gray(image_a);
[g,gof] = fit(z,image_b,equ,'Lower',[0.7,((hght+1)/2)-(0.25*hght),1,0],'Upper',[1.3,((hght+1)/2)+(0.25*hght),hght*0.8,0.3],'StartPoint',[1,(hght+1)*0.5,hght*0.3,0]);

% pull out the coefficients
coeff = coeffvalues(g);

if gof.rsquare > gof_thresh
    % if the gof is high enough, convert the c coeff to FWHM
    FWHM = (coeff(3)/sqrt(2)) * sqrt(8*log(2));
else
    % if it doesn't, log a dummy variable that will be excluded
    FWHM = -1;
end