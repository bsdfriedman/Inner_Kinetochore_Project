function [output] = IKK_sim_convert_ver2(infile,beads_inward,every_n_timesteps)

% The beads_inward input denotes the number of beads inward from the
% centromere you want to be labeled. For example, if it is 1, the two beads
% adjacent to each centromere will be labeled white. You can create
% multiple color files by making a larger array (E.g. [2 5] will make two
% color files with bead 2 labeled in one, and bead 5 labeled in the other).

%% Initialize
% log the number of masses, springs, and hinges that will be used in the calculation
[~,mass_num] = system(sprintf('grep "mass " %s | wc -l',infile));
mass_num = str2double(mass_num);
[~,spring_num] = system(sprintf('grep "spring " %s | wc -l',infile));
spring_num = str2double(spring_num);
[~,hinge_num] = system(sprintf('grep "hinge " %s | wc -l',infile));
hinge_num = str2double(hinge_num);

% log the time steps that will be used in the calculation
[~,times] = system(sprintf('grep "Time " %s | gawk ''!(NR%%%d)''',infile,every_n_timesteps));
time_name_temp = strsplit(times,'\n')';
time_name_temp = time_name_temp(~cellfun('isempty',time_name_temp));

% crop out unnecessary time points
if mod(size(time_name_temp,1),2) == 0
    time_name = time_name_temp(((size(time_name_temp,1)/2)+1):size(time_name_temp,1));
else
    time_name = time_name_temp(((size(time_name_temp,1)+1)/2):size(time_name_temp,1));
end

% parameters
mass_mass = 3.38889e-020; % mass of a standard bead
mass_multiplier = 10^10; % pinning mass
ignore = 0; % how many time points you want to ignore
mt_radius = 1.25*10^-7; % radius of microtubules
pin_2_cen = 1.5*10^-7; %distance from pinned beads to the centromere
kk_base = 8*10^-7; %distance between kinetochores
spring_rest = 1e-008; % spring distance at rest
spring_const = 0.226195; % standard spring constant
hinge_const = 4.0715e-012; % standard hinge constant

% preallocate
cen_idx = zeros([16 1]);
MT_idx = zeros([16 1]);
coords = zeros([mass_num 3 size(time_name,1)]);
coords_t0 = zeros([mass_num 3]);
bead_prop = zeros([mass_num 2]); % bead mass and colors for the reprinting (commented out)
springs = zeros([spring_num 2]);
hinges = zeros([hinge_num 3]);

% counters
cen_counter = 1;
time_counter = 1;
MT_counter = 1;
prop_counter = 1;
spring_counter = 1;
hinge_counter = 1;

%% Create the timepoint files, log colors, and log coordinates
% open up the file
fid_in = fopen(infile);

% assign tline so that the lines can be looped through
tline = fgetl(fid_in);

while ischar(tline) == 1
    if size(strfind(tline,'mass '),1) ~= 0
        % split the string into pieces to parse coordinates
        b = strsplit(tline);
        bead_prop(prop_counter,1) = str2double(b{4});
        if size(b,2) == 7
            % rows without an eighth entry correpond to a color of 1 (red)
            mass_color = 1;
            bead_prop(prop_counter,2) = 1;
        elseif size(b,2) == 8
            % rows with an eighth entry have their color logged appropriately
            mass_color = str2double(b{8});
            bead_prop(prop_counter,2) = str2double(b{8});
        else
        end
        if mass_color == 2 && str2double(b{4}) < mass_mass*mass_multiplier
            % log the appropriate cen index numbers
            cen_idx(cen_counter,1) = str2double(b{3});
            % increase the counter by 1
            cen_counter = cen_counter+1;
        elseif mass_color == 3
            % log the appropriate MT index numbers
            MT_idx(MT_counter,1) = str2double(b{3});
            % increase the counter by 1
            MT_counter = MT_counter+1;
        else
        end
        coords_t0(prop_counter,1:3) = [str2double(b{5}) str2double(b{6}) str2double(b{7})];
        % increase the counter by 1
        prop_counter = prop_counter + 1;
        % loop to the next line
        tline = fgetl(fid_in);
    elseif size(strfind(tline,'spring '),1) ~= 0
        % split the string into pieces to parse spring connections
        c = strsplit(tline);
        % log the spring connections
        springs(spring_counter,1:2) = [str2double(c{3}) str2double(c{4})];
        % increase the counter by 1
        spring_counter = spring_counter + 1;
        % loop to the next line
        tline = fgetl(fid_in);
    elseif size(strfind(tline,'hinge '),1) ~= 0
        % split the string into pieces to parse hinge connections
        d = strsplit(tline);
        % log the spring connections
        hinges(hinge_counter,1:3) = [str2double(d{3}) str2double(d{4}) str2double(d{5})];
        % increase the counter by 1
        hinge_counter = hinge_counter + 1;
        % loop to the next line
        tline = fgetl(fid_in);
    elseif size(strfind(tline,'MassColors'),1) ~= 0
        mass_color_list = zeros([mass_num 1]);
        if beads_inward > cen_idx(1,1)
            % if the chosen position exceeds the maximum possible value
            error('This position is too far from the centromere')
        else
        end
        for z = 1:mass_num
            % loop to the next line
            tline = fgetl(fid_in);
            % log the color
            mass_color_list(z,1) = str2double(tline);
        end
        % loop to the next line
        tline = fgetl(fid_in);
    elseif max(strcmp(tline,time_name)) == 1
        % loop to the next line
        tline = fgetl(fid_in);
        for z = 1:mass_num
            coords(z,1:3,time_counter) = str2double(strsplit(tline));
            % loop to the next line
            tline = fgetl(fid_in);
        end
        % increase the time counter by 1
        time_counter = time_counter + 1;
    else
        % loop to the next line
        tline = fgetl(fid_in);
    end
end

% close all the files
fclose('all');

%% Creating the timepoints file
% assign the name of the new file that is to be created
file_root_name = strsplit(infile,'.');
filename = sprintf('%s_timepoints.txt',file_root_name{1});
fid_out = fopen(filename,'w');

for z = 1:size(time_name,1)
    % pring the time line
    fprintf(fid_out,'%s\r\n',time_name{z,1});
    % print the first kinetochore
    for h = 1:size(coords,1)
        fprintf(fid_out,'%0.6f %0.6f %0.6f\r\n',(coords(h,1,z)+pin_2_cen)*10^6,coords(h,2,z)*10^6,coords(h,3,z)*10^6);
    end
    % print the second kinetochore
    for h = 1:size(coords,1)
        fprintf(fid_out,'%0.6f %0.6f %0.6f\r\n',(-(coords(h,1,z)+pin_2_cen)-kk_base)*10^6,coords(h,2,z)*10^6,coords(h,3,z)*10^6);
    end
    % print a blank line
    if z < size(time_name,1)
        fprintf(fid_out,'\r\n');
    end
end

% close all the files
fclose('all');

%% Creating the color files
% preallocate the table for the recolored beads
recolor_idx = zeros([32 1 ]);

% choose the beads that are to be recolored
for q = 1:size(beads_inward,2)
    for z = 1:size(cen_idx,1)
        recolor_idx(z*2-1,q) = cen_idx(z,1)-beads_inward(q);
        recolor_idx(z*2,q) = cen_idx(z,1)+beads_inward(q);
    end
end

% preallocate
mass_color_list_new = zeros([size(mass_color_list,1) 1 size(recolor_idx,2)]);

for q = 1:size(recolor_idx,2)
    % create the colors file
    filename = sprintf('%s_colors_%d.txt',file_root_name{1},beads_inward(q));
    fid_out = fopen(filename,'w');
    
    % assign the mass colors
    mass_color_list_new(:,:,q) = mass_color_list;
    
    % change the appropriate colors in the mass color list
    for z = 1:size(recolor_idx,1)
        mass_color_list_new(mass_num-recolor_idx(z,q),1,q) = 4;
    end
    
    for y = 1:2
        for z = 1:size(mass_color_list_new,1)
            % print out the color names into the color file
            fprintf(fid_out,'%d\r\n',mass_color_list_new(z,1,q));
        end
    end
    
    % close all the files
    fclose('all');
end

% make a color file for the MT beads
filename = sprintf('%s_colors_MT.txt',file_root_name{1});
fid_out = fopen(filename,'w');

% assign the mass colors
mass_color_list_MT = mass_color_list;

% change the appropriate colors in the mass color list
for z = 1:size(MT_idx,1)
    mass_color_list_MT(mass_num-MT_idx(z,1)) = 4;
end

for y = 1:2
    for z = 1:size(mass_color_list_MT,1)
        % print out the color names into the color file
        fprintf(fid_out,'%d\r\n',mass_color_list_MT(z,1));
    end
end

% close all the files
fclose('all');

% THIS CODE IS USED TO TEST THE PROGRAM AND HAS BEEN COMMENTED OUT
for q = 1:size(beads_inward,2)+1
    
    if q == size(beads_inward,2)+1
        % name the MT file if this is the last loop through
        filename = sprintf('%s_reprint_MT.out',file_root_name{1});
        fid_out = fopen(filename,'w');
    else
        % name the file
        filename = sprintf('%s_reprint_%d.out',file_root_name{1},beads_inward(q));
        fid_out = fopen(filename,'w');
    end
    
    % print out the intro
    fprintf(fid_out,'meta temperature_Celsius 25\r\nmeta viscosity_centiPoise 1\r\nmeta effective_damping_radius 8e-009\r\nmeta dna_modulus_gigaPascal 2\r\nmeta dna_radius_nanometers 0.6\r\nmeta damping_radius_factor 0.8\r\nstructure {\r\n  random_force 2.78554e-011\r\n  mass_damping 4.44973e+009\r\n  mass_radius 4.5e-009\r\n  time_step 2e-009\r\n  collision_spring_constant 0.0565487\r\n  spring_damping_factor 0\r\n  random_number_seed 42\r\n  color 1\r\n');
    % reprint the masses
    for z = 1:mass_num
        fprintf(fid_out,'  mass %d\t%.6g\t%.6g %.6g %.6g %d\r\n',z-1,bead_prop(z,1),(coords_t0(z,1)+pin_2_cen),coords_t0(z,2),coords_t0(z,3),bead_prop(z,2));
    end
    % print the duplicate masses
    for z = 1:mass_num
        fprintf(fid_out,'  mass %d\t%.6g\t%.6g %.6g %.6g %d\r\n',z+mass_num-1,bead_prop(z,1),(-(coords_t0(z,1)+pin_2_cen)-kk_base),coords_t0(z,2),coords_t0(z,3),bead_prop(z,2));
    end
    % reprint the springs
    for z = 1:spring_num
        fprintf(fid_out,'  spring %d %d %.1g %.6g\r\n',springs(z,1),springs(z,2),spring_rest,spring_const);
    end
    % print the duplicate springs
    for z = 1:spring_num
        fprintf(fid_out,'  spring %d %d %.1g %.6g\r\n',springs(z,1)+mass_num,springs(z,2)+mass_num,spring_rest,spring_const);
    end
    % reprint the hinges
    for z = 1:hinge_num
        fprintf(fid_out,'  hinge %d %d %d %.6g\r\n',hinges(z,1),hinges(z,2),hinges(z,3),hinge_const);
    end
    % print the duplicate hinges
    for z = 1:hinge_num
        fprintf(fid_out,'  hinge %d %d %d %.6g\r\n',hinges(z,1)+mass_num,hinges(z,2)+mass_num,hinges(z,3)+mass_num,hinge_const);
    end
    % print out the end of the cfg structure, then move down one line and print the MassColors line
    fprintf(fid_out,'}\r\n\r\nMassColors\r\n');
    
    for y = 1:2
        for z = 1:size(mass_color_list,1)
            if q == size(beads_inward,2)+1
                % print out the color names into the file
                fprintf(fid_out,'%d\r\n',mass_color_list_MT(z,1));
            else
                % print out the color names into the file
                fprintf(fid_out,'%d\r\n',mass_color_list_new(z,1,q));
            end
        end
    end
    
    % move down one line
    fprintf(fid_out,'\r\n');
    
    % print out the times
    for z = 1:size(time_name,1)
        fprintf(fid_out,'%s\r\n',time_name{z,1});
        % reprint out all the coordinates
        for h = 1:mass_num
            fprintf(fid_out,'%0.6e %0.6e %0.6e\r\n',(coords(h,1,z)+pin_2_cen),coords(h,2,z),coords(h,3,z));
        end
        % print the duplicate coordinates
        for h = 1:mass_num
            fprintf(fid_out,'%0.6e %0.6e %0.6e\r\n',(-(coords(h,1,z)+pin_2_cen)-kk_base),coords(h,2,z),coords(h,3,z));
        end
        if z < size(time_name,1)
            % move down one line if it is not the last time point
            fprintf(fid_out,'\r\n');
        else
        end
    end
end

% close all the files
fclose('all');

%% Performing the simulation measurements
% flip the coordinates so that they are right side up
coords = flipud(coords);

% preallocate the table for the recolored beads
bead_meas_idx = zeros([32 1]);

% choose the beads that are to be recolored
for q = 1:size(beads_inward,2)
    for z = 1:size(cen_idx,1)
        bead_meas_idx(z*2-1,q) = cen_idx(z,1)-beads_inward(q);
        bead_meas_idx(z*2,q) = cen_idx(z,1)+beads_inward(q);
    end
end

% preallocate the coordinate tables
cen_coords = zeros([16 3 size(beads_inward,2) size(time_name,1)]);
bead_coords = zeros([32 3 size(beads_inward,2) size(time_name,1)]);
MT_coords = zeros([16 3 size(time_name,1)]);

% log the coordinates for all beads
for y = 1:size(time_name,1)
    for q = 1:size(beads_inward,2)
        for z = 1:16
            cen_coords(z,1:3,q,y) = coords(cen_idx(z,1)+1,1:3,y);
            bead_coords(z*2-1,1:3,q,y) = coords(bead_meas_idx(z*2-1,q)+1,1:3,y);
            bead_coords(z*2,1:3,q,y) = coords(bead_meas_idx(z*2,q)+1,1:3,y);
        end
    end
end

% log the coordinates for MT beads
for y = 1:size(time_name,1)
    for z = 1:16
        MT_coords(z,1:3,y) = coords(MT_idx(z,1)+1,1:3,y);
    end
end

% preallocate the measurements
cen_bead_dist_table = zeros([32 1 size(beads_inward,2) size(time_name,1)]);
cen_bead_xpos_table = zeros([32 1 size(beads_inward,2) size(time_name,1)]);

% get the measurements
for y = 1:size(time_name,1)
    for q = 1:size(beads_inward,2)
        for z = 1:16
            % calculate the 3D distances between beads
            cen_bead_dist_table(z*2-1,1,q,y) = norm(cen_coords(z,1:3,q,y)-bead_coords(z*2-1,1:3,q,y));
            cen_bead_dist_table(z*2,1,q,y) = norm(cen_coords(z,1:3,q,y)-bead_coords(z*2,1:3,q,y));
            
            % calculate the position along the axis
            cen_bead_xpos_table(z*2-1,1,q,y) = cen_coords(z,1,q,y)-bead_coords(z*2-1,1,q,y);
            cen_bead_xpos_table(z*2,1,q,y) = cen_coords(z,1,q,y)-bead_coords(z*2,1,q,y);
        end
    end
end

% set up the output distances
cen_bead_dist = cell([size(beads_inward,2)+1 2]);
cen_bead_dist{1,1} = 'Bead Position';
cen_bead_dist{1,2} = 'Distance(nm)';

% set up the output x-positions
cen_bead_xpos = cell([size(beads_inward,2)+1 2]);
cen_bead_xpos{1,1} = 'Bead Position';
cen_bead_xpos{1,2} = 'X-Position(nm)';

if 1+ignore > size(cen_bead_dist_table,4) || 1+ignore > size(cen_bead_xpos_table,4)
    error('You cannot ignore that many timepoints')
end

% convert everything to nm
cen_bead_dist_table = cen_bead_dist_table*10^9;
cen_bead_xpos_table = cen_bead_xpos_table*10^9;

% calculate the mean values over all time points
for z = 1:size(beads_inward,2)
    cen_bead_dist{z+1,2} = mean(mean(mean(cen_bead_dist_table(:,:,z,1+ignore:size(cen_bead_dist_table,4)))));
    cen_bead_xpos{z+1,2} = mean(mean(mean(cen_bead_xpos_table(:,:,z,1+ignore:size(cen_bead_xpos_table,4)))));
    cen_bead_dist{z+1,1} = beads_inward(z);
    cen_bead_xpos{z+1,1} = beads_inward(z);
end

% set up the output x-positions
xpos_percent = cell([size(beads_inward,2)+1 4]);
xpos_percent{1,1} = 'Bead Position';
xpos_percent{1,2} = 'Percent DNA';
xpos_percent{1,3} = 'Percent SPB';
xpos_percent{1,4} = 'Percent Same';

% calculate the position percentages
for z = 1:size(beads_inward,2)
    cen_measure = cen_bead_xpos_table(:,:,z,1+ignore:size(cen_bead_xpos_table,4));
    xpos_percent{z+1,2} = sum(sum(cen_measure > 0))/(size(cen_measure,1)*size(cen_measure,2)*size(cen_measure,3)*size(cen_measure,4));
    xpos_percent{z+1,3} = sum(sum(cen_measure < 0))/(size(cen_measure,1)*size(cen_measure,2)*size(cen_measure,3)*size(cen_measure,4));
    xpos_percent{z+1,4} = sum(sum(cen_measure == 0))/(size(cen_measure,1)*size(cen_measure,2)*size(cen_measure,3)*size(cen_measure,4));
    xpos_percent{z+1,1} = beads_inward(z);
end

% preallocate the measurements
MT_cen_xpos_table = zeros([16 1 size(time_name,1)]);

% get the measurements
for y = 1:size(time_name,1)
    for z = 1:16
        % calculate the position along the axis
        MT_cen_xpos_table(z,1,y) = MT_coords(z,1,y)-cen_coords(z,1,1,y);
    end
end

% convert to nm
MT_cen_xpos_table = MT_cen_xpos_table*10^9;

% set up the output x-positions
MTpos_percent = cell([2 3]);
MTpos_percent{1,1} = 'Percent DNA';
MTpos_percent{1,2} = 'Percent SPB';
MTpos_percent{1,3} = 'Percent Same';

% calculate the position percentages
MT_measure = MT_cen_xpos_table(:,:,1+ignore:size(MT_cen_xpos_table,3));
MTpos_percent{2,1} = sum(sum(MT_measure > 0))/(size(MT_measure,1)*size(MT_measure,2)*size(MT_measure,3));
MTpos_percent{2,2} = sum(sum(MT_measure < 0))/(size(MT_measure,1)*size(MT_measure,2)*size(MT_measure,3));
MTpos_percent{2,3} = sum(sum(MT_measure == 0))/(size(MT_measure,1)*size(MT_measure,2)*size(MT_measure,3));

% preallocate the bead radius table
bead_rad_table = zeros([32 size(beads_inward,2) size(time_name,1)]);

% calculate the bead radii
for q = 1:size(time_name,1)
    for z = 1:size(beads_inward,2)
        for h = 1:16
            bead_rad_table(h*2-1,z,q) = norm(bead_coords(h*2-1,2:3,z,q));
            bead_rad_table(h*2,z,q) = norm(bead_coords(h*2,2:3,z,q));
        end
    end
end

% set up the output for bead radius
bead_radius = cell([size(beads_inward,2)+2 2]);
bead_radius{1,1} = 'Bead Position';
bead_radius{1,2} = 'Average Radius (nm)';
bead_radius{2,1} = 'MT';
bead_radius{2,2} = mt_radius*10^9;

% find the average radii for each chosen bead
for z = 1:size(beads_inward,2)
    bead_radius{z+2,1} = beads_inward(z);
    bead_radius{z+2,2} = mean(mean(bead_rad_table(:,z,:)))*10^9;
end

% create the outputs
output.cen_bead_dist_table = cen_bead_dist_table;
output.cen_bead_xpos_table = cen_bead_xpos_table;
output.MT_cen_xpos_table = MT_cen_xpos_table;
output.bead_rad_table = bead_rad_table;
output.cen_bead_dist = cen_bead_dist;
output.cen_bead_xpos = cen_bead_xpos;
output.xpos_percent = xpos_percent;
output.MTpos_percent = MTpos_percent;
output.bead_radius = bead_radius;

