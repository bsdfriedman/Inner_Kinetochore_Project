function ndc80_matching_theta

% this code builds two rings where the corresponding beads are always 70 nm
% apart, and the second ring, representing ame1 gets larger in diameter and
% closer to the smaller ring, representing ndc80, as the kinetochore
% swivels up by a chosen angle, in each iteration

% upon revision, it now makes two kinetochores

%% Get the bead coordinates for the circles

% input
pivot_angle = 5; % the angle in degrees that the large ring pivots upwards

% parameters
ndc80_radius_nm = 130; % N-ndc80 radius in nm
theta_1 = (2*pi)/16; % angle of rotation to make circles of 16 beads
kin_length = 30; % length of the kinetochore in nm
num_kin = 16; % number of kinetochores
create_outfile = 1; % this is 0 if you don't want an out file and 1 if you do
mass_mass = 3.38889e-020; % mass of a standard bead
mass_multiplier = 10^10; % pinning mass
kin_kin_dist = 10^-6; %distance between kinetochores in m

% select all the angles for the pivot
theta_2_deg = (0:(90/pivot_angle))*pivot_angle;
theta_2 = deg2rad(theta_2_deg);% the table of pivot angles

% calculate the x displacement and radii
x_disp_nm = cos(theta_2)*kin_length;
rad_nm = sin(theta_2)*kin_length + ndc80_radius_nm;

% add the ndc80 ring onto the two tables
x_disp_nm = cat(2,0,x_disp_nm);
rad_nm = cat(2,ndc80_radius_nm,rad_nm);

% convert the values from nm to meters
x_disp_m = x_disp_nm*10^-9;
rad_m = rad_nm*10^-9;

% preallocate
bead_coords_m = zeros([num_kin*size(rad_nm,2) 3]);

% assign the y and z coordinates in m, for the chromoShake out file
for y = 1:size(rad_nm,2)
    for z = 1:num_kin
        bead_coords_m(z+(num_kin*(y-1)),1:3) = [1 0 0; 0 cos(theta_1*(z-1)) -sin(theta_1*(z-1)); 0 sin(theta_1*(z-1)) cos(theta_1*(z-1))]*[0; rad_m(y); 0];
    end
end

% assign the x_coordinates
for z = 1:size(x_disp_nm,2)
    bead_coords_m(((z-1)*num_kin+1):z*num_kin,1) = x_disp_m(z);
end

% flip the coordinates file to be input into the out/timepoints file
bead_coords_m = flipud(bead_coords_m);

%% Create the timepoint file
% create the time file
timefile_name = sprintf('ndc80_theta_%d_%d_%d_timepoints.txt',min(theta_2_deg),max(theta_2_deg),theta_2_deg(2)-theta_2_deg(1));
fid_out = fopen(timefile_name,'w');

% print the time line
fprintf(fid_out,'Time 0\r\n');
% print out the coordinates for each bead
for z = 1:size(bead_coords_m)
    fprintf(fid_out,'%0.6f %0.6f %0.6f\r\n',bead_coords_m(z,1)*10^6,bead_coords_m(z,2)*10^6,bead_coords_m(z,3)*10^6);
end
% print out the coordinates for the second set of beads
for z = 1:size(bead_coords_m)
    fprintf(fid_out,'%0.6f %0.6f %0.6f\r\n',(-bead_coords_m(z,1)+kin_kin_dist)*10^6,bead_coords_m(z,2)*10^6,bead_coords_m(z,3)*10^6);
end

% close the file
fclose('all');

%% Create the color files
% preallocate
colors_new = zeros([size(bead_coords_m,1) size(rad_m,2)]);

for z = 1:size(rad_m,2)
    
    if z == 1
        % create the time file
        colorfile_name = sprintf('ndc80_theta_%d_%d_%d_colors_rndc80.txt',min(theta_2_deg),max(theta_2_deg),theta_2_deg(2)-theta_2_deg(1));
        fid_out = fopen(colorfile_name,'w');
    else
        % create the time file
        colorfile_name = sprintf('ndc80_theta_%d_%d_%d_colors_r%d.txt',min(theta_2_deg),max(theta_2_deg),theta_2_deg(2)-theta_2_deg(1),theta_2_deg(z-1));
        fid_out = fopen(colorfile_name,'w');
    end
    
    % assign the base colors
    colors_new(:,:,z) = 5; % this makes the normal beads white
    colors_new(((z-1)*num_kin+1):z*num_kin,z) = 4; % makes this ring fluorescent
    
    % move the colors to a new table to print them
    colors_log = colors_new(:,z);
    
    % flip the color file to match the timepoints
    colors_log = flipud(colors_log);
    
    % print out the colors
    for h = 1:2
        for y = 1:size(colors_log)
            fprintf(fid_out,'%d\r\n',colors_log(y,1));
        end
    end
    
    % close the file
    fclose('all');
end

%% Create the out file if we want it
if create_outfile == 1
    for z = 1:size(rad_nm,2)
        if z == 1
            % create the out file
            outfile_name = sprintf('ndc80_theta_%d_%d_%d_rndc80.out',min(theta_2_deg),max(theta_2_deg),theta_2_deg(2)-theta_2_deg(1));
            fid_out = fopen(outfile_name,'w');
        else
            % create the out file
            outfile_name = sprintf('ndc80_theta_%d_%d_%d_r%d.out',min(theta_2_deg),max(theta_2_deg),theta_2_deg(2)-theta_2_deg(1),theta_2_deg(z-1));
            fid_out = fopen(outfile_name,'w');
        end
        
        % print the intro
        fprintf(fid_out,'meta temperature_Celsius 25\r\nmeta viscosity_centiPoise 1\r\nmeta effective_damping_radius 8e-009\r\nmeta dna_modulus_gigaPascal 2\r\nmeta dna_radius_nanometers 0.6\r\nmeta damping_radius_factor 0.8\r\nstructure {\r\n  random_force 2.78554e-011\r\n  mass_damping 4.44973e+009\r\n  mass_radius 4.5e-009\r\n  time_step 2e-009\r\n  collision_spring_constant 0.0565487\r\n  spring_damping_factor 0\r\n  random_number_seed 42\r\n  color 1\r\n');
        
        % temporarily flip the coordinates and color files to create the masses
        bead_coords_m = flipud(bead_coords_m);
        
        % print the masses
        for y = 1:size(bead_coords_m,1)
            fprintf(fid_out,'  mass %d\t%.6g\t%.6g %.6g %.6g %d\r\n',y-1,mass_mass*mass_multiplier,bead_coords_m(y,1),bead_coords_m(y,2),bead_coords_m(y,3),colors_new(y,z));
        end
        % print the duplicate masses
        for y = 1:size(bead_coords_m,1)
            fprintf(fid_out,'  mass %d\t%.6g\t%.6g %.6g %.6g %d\r\n',y-1+size(bead_coords_m,1),mass_mass*mass_multiplier,-bead_coords_m(y,1)+kin_kin_dist,bead_coords_m(y,2),bead_coords_m(y,3),colors_new(y,z));
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
        for h = 1:2
            for y = 1:size(colors_log)
                fprintf(fid_out,'%d\r\n',colors_log(y,1));
            end
        end
        
        % move down one line and print the 'Time 0' line
        fprintf(fid_out,'\r\nTime 0\r\n');
        
        % print the coordinates for the timepoint
        for y = 1:size(bead_coords_m,1)
            fprintf(fid_out,'%0.6e %0.6e %0.6e\r\n',bead_coords_m(y,1),bead_coords_m(y,2),bead_coords_m(y,3));
        end
        % print the duplicate coordinates for the timepoint
        for y = 1:size(bead_coords_m,1)
            fprintf(fid_out,'%0.6e %0.6e %0.6e\r\n',-bead_coords_m(y,1)+kin_kin_dist,bead_coords_m(y,2),bead_coords_m(y,3));
        end
        
        % close the file
        fclose('all');
    end
else
end