function IKK_sim_convert_horiz_cen(infile,every_n_timesteps)

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

filename = sprintf('%s_colors_horiz_cen.txt',file_root_name{1});
fid_out = fopen(filename,'w');

% assign the mass colors
mass_color_list_new = mass_color_list;

% preallocate the table for the recolored beads
recolor_idx = zeros([48 1]);

% mark the beads to be recolored
for z = 1:size(cen_idx,1)
    recolor_idx(z*3-2,1) = cen_idx(z,1)-1;
    recolor_idx(z*3-1,1) = cen_idx(z,1);
    recolor_idx(z*3,1) = cen_idx(z,1)+1;
end

% change the appropriate colors in the mass color list
for z = 1:size(recolor_idx,1)
    mass_color_list_new(mass_num-recolor_idx(z),1) = 4;
end

for y = 1:2
    for z = 1:size(mass_color_list_new,1)
        % print out the color names into the color file
        fprintf(fid_out,'%d\r\n',mass_color_list_new(z,1));
    end
end

% close all the files
fclose('all');

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

% name the file
filename = sprintf('%s_reprint_horiz_cen.out',file_root_name{1});
fid_out = fopen(filename,'w');

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
        % print out the color names into the file
        fprintf(fid_out,'%d\r\n',mass_color_list_new(z,1));
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

% close all the files
fclose('all');

