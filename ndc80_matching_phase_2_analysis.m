function [p_vals] = ndc80_matching_phase_2_analysis(cells_FWHM,sim_FWHM)

% parameters
min_angle = 0; % minimum radius in nm
max_angle = 90; % minimum radius in nm
angle_step_size = 5; % increase in radius with each step

% load the step-sizes to populaiton the table
rad_nm = (0:((max_angle-min_angle)/angle_step_size))*angle_step_size+min_angle;

% preallocate
p_vals = cell([size(sim_FWHM,1)+1 2]);
p_vals(1,1:2) = {'Circle Radius','kstest2 p_value'};
for z = 1:size(sim_FWHM,1)-1
    p_vals{z+1,1} = rad_nm(z);
end
p_vals(size(sim_FWHM,1)+1,1) = {'ndc80'};

% loop through and perform a kstest2 on the data
for z = 1:size(sim_FWHM,1)
    % remove all negative values
    sim_data = sim_FWHM(z,:);
    sim_data = sim_data(sim_data>0);
    
    % perform a kstest2 and pull out the p value
    [~,p] = kstest2(cells_FWHM,sim_data);
    
    % store the p_value in the table
    p_vals{z+1,2} = p;
end

