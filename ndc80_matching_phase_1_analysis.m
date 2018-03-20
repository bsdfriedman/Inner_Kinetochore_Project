function [p_vals] = ndc80_matching_phase_1_analysis(cells_FWHM,sim_FWHM)

% parameters
min_radius = 100; % minimum radius in nm
max_radius = 140; % minimum radius in nm
radius_step_size = 1; % increase in radius with each step

% load the step-sizes to populaiton the table
rad_nm = (0:((max_radius-min_radius)/radius_step_size))*radius_step_size+min_radius;

% preallocate
p_vals = cell([size(sim_FWHM,1)+1 2]);
p_vals(1,1:2) = {'Circle Radius','kstest2 p_value'};
for z = 1:size(sim_FWHM,1)
    p_vals{z+1,1} = rad_nm(z);
end

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

