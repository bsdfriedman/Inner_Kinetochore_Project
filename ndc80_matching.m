function ndc80_matching

% this code creates a series of concentric circles to see which one matches
% the ndc80 distribution the best

%% Get the bead coordinates for the circles
% parameters
num_kin = 16; % number of kinetochores
theta = (2*pi)/16; % angle of rotation
min_radius = 100; % minimum radius in nm
max_radius = 140; % minimum radius in nm
radius_step_size = 1; % increase in radius with each step
create_outfile = 0; % this is 0 if you don't want an out file and 1 if you do
mass_mass = 3.38889e-020; % mass of a standard bead
mass_multiplier = 10^10; % pinning mass

% choose the radii of the circles in nm
rad_nm = (0:((max_radius-min_radius)/radius_step_size))*radius_step_size+min_radius;
% convert the radii into meters
rad_m = rad_nm*10^-9;

% preallocate
bead_coords_m = zeros([num_kin*size(rad_nm,2) 3]);

% assign the coordinates in m, for the chromoShake out file
for y = 1:size(rad_nm,2)
    for z = 1:num_kin
        bead_coords_m(z+(num_kin*(y-1)),1:3) = [1 0 0; 0 cos(theta*(z-1)) -sin(theta*(z-1)); 0 sin(theta*(z-1)) cos(theta*(z-1))]*[0; rad_m(y); 0];
    end
end

% flip the coordinates file to be input into the out/timepoints file
bead_coords_m = flipud(bead_coords_m);

%% Create the timepoint file
% create the time file
timefile_name = sprintf('ndc80_circles_%d_%d_%d_timepoints.txt',min(rad_nm),max(rad_nm),rad_nm(2)-rad_nm(1));
fid_out = fopen(timefile_name,'w');

% print the time line
fprintf(fid_out,'Time 0\r\n');
% print out the coordinates for each bead
for z = 1:size(bead_coords_m)
    fprintf(fid_out,'%0.6f %0.6f %0.6f\r\n',bead_coords_m(z,1)*10^6,bead_coords_m(z,2)*10^6,bead_coords_m(z,3)*10^6);
end

% close the file
fclose('all');

%% Create the color files
% preallocate
colors_new = zeros([size(bead_coords_m,1) size(rad_m,2)]);

for z = 1:size(rad_m,2)
    % create the time file
    colorfile_name = sprintf('ndc80_circles_%d_%d_%d_colors_r%d.txt',min(rad_nm),max(rad_nm),rad_nm(2)-rad_nm(1),rad_nm(z));
    fid_out = fopen(colorfile_name,'w');
    
    % assign the base colors
    colors_new(:,:,z) = 5; % this makes the normal beads white
    colors_new(((z-1)*num_kin+1):z*num_kin,z) = 4; % makes this ring fluorescent
    
    % move the colors to a new table to print them
    colors_log = colors_new(:,z);
    
    % flip the color file to match the timepoints
    colors_log = flipud(colors_log);
    
    % print out the colors
    for y = 1:size(colors_log)
        fprintf(fid_out,'%d\r\n',colors_log(y,1));
    end
    
    % close the file
    fclose('all');
end

%% Create the out file if we want it
if create_outfile == 1
    for z = 1:size(rad_nm,2)
        % create the out file
        outfile_name = sprintf('ndc80_circles_%d_%d_%d_r%d.out',min(rad_nm),max(rad_nm),rad_nm(2)-rad_nm(1),rad_nm(z));
        fid_out = fopen(outfile_name,'w');
        
        % print the intro
        fprintf(fid_out,'meta temperature_Celsius 25\r\nmeta viscosity_centiPoise 1\r\nmeta effective_damping_radius 8e-009\r\nmeta dna_modulus_gigaPascal 2\r\nmeta dna_radius_nanometers 0.6\r\nmeta damping_radius_factor 0.8\r\nstructure {\r\n  random_force 2.78554e-011\r\n  mass_damping 4.44973e+009\r\n  mass_radius 4.5e-009\r\n  time_step 2e-009\r\n  collision_spring_constant 0.0565487\r\n  spring_damping_factor 0\r\n  random_number_seed 42\r\n  color 1\r\n');
        
        % temporarily flip the coordinates and color files to create the masses
        bead_coords_m = flipud(bead_coords_m);
        
        % print the masses
        for y = 1:size(bead_coords_m,1)
            fprintf(fid_out,'  mass %d\t%.6g\t%.6g %.6g %.6g %d\r\n',y-1,mass_mass*mass_multiplier,bead_coords_m(y,1),bead_coords_m(y,2),bead_coords_m(y,3),colors_new(y,z));
        end
        
        % flip the coordinates file back to create the timepoint
        bead_coords_m = flipud(bead_coords_m);
        
        % print out the end of the cfg structure, then move down one line and print the MassColors line
        fprintf(fid_out,'}\r\n\r\nMassColors\r\n');
        
        % move the colors to a new table to print them
        colors_log = colors_new(:,z);
        
        % flip the color file to match the timepoints
        colors_log = flipud(colors_log);
        
        % print out the colors
        for y = 1:size(colors_log)
            fprintf(fid_out,'%d\r\n',colors_log(y,1));
        end
        
        % move down one line and print the 'Time 0' line
        fprintf(fid_out,'\r\nTime 0\r\n');
        
        % print the coordinates for the timepoint
        for y = 1:size(bead_coords_m,1)
            fprintf(fid_out,'%0.6e %0.6e %0.6e\r\n',bead_coords_m(y,1),bead_coords_m(y,2),bead_coords_m(y,3));
        end
        
        % close the file
        fclose('all');
    end
else
end
