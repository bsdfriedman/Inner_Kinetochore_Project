function [output] = IKK_chromoShake_bead_position(infile, num_timesteps,beads_inward)

% log the time steps that will be used in the calculation
[~,mass_num] = system(sprintf('grep "mass " %s | wc -l',infile));
mass_num = str2double(mass_num);

% log the time steps that will be used in the calculation
[~,times] = system(sprintf('grep "Time " %s | tail -%d',infile,num_timesteps));
time_name = strsplit(times,'\n')';
time_name = time_name(~cellfun('isempty',time_name));

% parameters
mass_mass = 3.38889e-020; % mass of a standard bead
mass_multiplier = 10^10; % pinning mass

% preallocate
cen_idx = zeros([16 1]);
coords = zeros([mass_num 3 size(time_name,1)]);

% counters
cen_counter = 1;
time_counter = 1;

% open up the file
fid_in = fopen(infile);

% assign tline so that the lines can be looped through
tline = fgetl(fid_in);

while ischar(tline) == 1
    if size(strfind(tline,'mass '),1) ~= 0
        % split the string into pieces to parse coordinates
        b = strsplit(tline);
        if size(b,2) == 7
            % rows without an eighth entry correpond to a color of 1 (red)
            mass_color = 1;
        elseif size(b,2) == 8
            % rows with an eighth entry have their color logged appropriately
            mass_color = str2double(b{8});
        else
        end
        if mass_color == 2 && str2double(b{4}) < mass_mass*mass_multiplier
            % log the appropriate cen index numbers
            cen_idx(cen_counter,1) = str2double(b{3})+1;
            % increase the counter by 1
            cen_counter = cen_counter+1;
        else
        end
        % loop to the next line
        tline = fgetl(fid_in);
    elseif max(strcmp(tline,time_name)) == 1
        % loop to the next line
        tline = fgetl(fid_in);
        for z = 1:mass_num
            % log the coordinates
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

% if the cen beads are chosen, remove duplicate values
bead_meas_idx = unique(bead_meas_idx,'rows');

% preallocate the coordinate tables
cen_coords = zeros([16 3 size(time_name,1)]);
bead_coords = zeros([32 3 size(time_name,1)]);

% log the coordinates for all beads
for q = 1:size(time_name,1)
    for z = 1:16
        cen_coords(z,1:3,q) = coords(cen_idx(z,1),1:3,q);
        bead_coords(z*2-1,1:3,q) = coords(bead_meas_idx(z*2-1,1),1:3,q);
        bead_coords(z*2,1:3,q) = coords(bead_meas_idx(z*2,1),1:3,q);
    end
end

% preallocate the measurements
cen_bead_dist_table = zeros([32 1 size(time_name,1)]);
cen_bead_xpos_table = zeros([32 1 size(time_name,1)]);

% get the measurements
for q = 1:size(time_name,1)
    for z = 1:16
        % calculate the 3D distances between beads
        cen_bead_dist_table(z*2-1,1,q) = norm(cen_coords(z,1:3,q)-bead_coords(z*2-1,1:3,q));
        cen_bead_dist_table(z*2,1,q) = norm(cen_coords(z,1:3,q)-bead_coords(z*2,1:3,q));
        
        % calculate the position along the axis
        cen_bead_xpos_table(z*2-1,1,q) = cen_coords(z,1,q)-bead_coords(z*2-1,1,q);
        cen_bead_xpos_table(z*2,1,q) = cen_coords(z,1,q)-bead_coords(z*2,1,q);
    end
end

% create the outputs
output.cen_bead_dist_table = cen_bead_dist_table;
output.cen_bead_xpos_table = cen_bead_xpos_table;
output.cen_bead_dist = mean(mean(cen_bead_dist_table));
output.cen_bead_xpos = mean(mean(cen_bead_xpos_table));

% calculate the percentage of times the bead is in a given position relative to the cen
output.xpos_percent_DNA = sum(sum(cen_bead_xpos_table > 0))/(size(time_name,1)*32);
output.xpos_percent_spindle = sum(sum(cen_bead_xpos_table < 0))/(size(time_name,1)*32);
output.xpos_percent_same = sum(sum(cen_bead_xpos_table == 0))/(size(time_name,1)*32);
