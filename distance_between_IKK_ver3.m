function [dist] = distance_between_IKK_ver3(a,b)

% this code caculates the xy distance between spots a and b
dist = sqrt((a(1)-b(1))^2 + (a(2)-b(2))^2);