function [plane] = plane_conv_IKK_ver3(input)
% Written by Brandon Friedman

% this code converts the plane number to a value between 1 and 7
plane = input - (floor(input/7)*7);
if plane == 0
    plane = 7;
else
end