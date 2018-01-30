function [theta] = angle_between_IKK_ver3(a, b)
% Written by Brandon Friedman
% this code calculates the angle between two spots relative to the positive x-axis

% get the difference in x and y
x = b(1) - a(1);
y = b(2) - a(2);

% use arctan to get the angle
theta = rad2deg(atan2(y,x));
end