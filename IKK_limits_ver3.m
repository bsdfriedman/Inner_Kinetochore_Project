function [limit, d_r1r2, d_g1g2] = IKK_limits_ver3(r1,r2,g1,g2,lim_r,lim_g,planesep,pixel_size)

% for limits, if any of the conditions are not met, the limits counter increases by 1
% the program will only perform calculations on the data if limit ends up being 0

% calculate the distance between the RFP spots
d_r1r2 = distance_between_IKK_ver3(r1,r2) * pixel_size;

% calculate the distance between the protein spots
d_g1g2 = distance_between_IKK_ver3(g1,g2) * pixel_size;

% make sure both distances are within range
a = d_r1r2 < lim_r(1);
b = d_r1r2 > lim_r(2);
c = d_g1g2 < lim_g(1);
d = d_g1g2 > lim_g(2);

% make sure none of the spots are duplicates
e = min(r1(1:3) == r2(1:3)) + min(g1(1:3) == g2(1:3));

% exclude based on tilt
f = abs(r1(3) - r2(3))> planesep;
g = abs(g1(3) - g2(3))> planesep;

limit = a+b+c+d+e+f+g;