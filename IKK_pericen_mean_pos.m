function [output] = IKK_pericen_mean_pos(infile,beads_inward,remove_first_n_timesteps)

% log the number of masses
[~,mass_num] = system(sprintf('grep "mass " %s | wc -l',infile));
mass_num = str2double(mass_num);
[~,time_num] = system(sprintf('grep "Time " %s | wc -l',infile));
time_num = str2double(time_num);

% log the time steps that will be excluded
[~,times] = system(sprintf('grep "Time " %s | head -%d',infile,remove_first_n_timesteps));
time_name = strsplit(times,'\n')';
time_name = time_name(1:remove_first_n_timesteps,1);

% parameters
mass_mass = 3.38889e-020; % mass of a standard bead
mass_multiplier = 10^10; % pinning mass

% preallocate
cen_idx = zeros([16 1]);
MT_idx = zeros([16 1]);
coords = zeros([mass_num 3 (time_num-size(time_name,1))]);

% counters
cen_counter = 1;
time_counter = 1;
MT_counter = 1;

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
    elseif size(strfind(tline,'Time '),1) ~= 0 && max(strcmp(tline,time_name)) == 0
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

% flip the coordinates rightside up
coords = flipud(coords);

% preallocate the table for the recolored beads
meas_idx = zeros([size(cen_idx,1)*size(beads_inward,2) 5]);

% choose the beads that are to be measured
for q = 1:size(beads_inward,2)
    for z = 1:size(cen_idx,1)
        meas_idx(z*2-1+((q-1)*size(cen_idx,1)),1) = -beads_inward(q);
        meas_idx(z*2-1+((q-1)*size(cen_idx,1)),2) = cen_idx(z,1)-beads_inward(q);
        meas_idx(z*2+((q-1)*size(cen_idx,1)),1) = beads_inward(q);
        meas_idx(z*2+((q-1)*size(cen_idx,1)),2) = cen_idx(z,1)+beads_inward(q);
    end
end

% remove all duplicate values
meas_idx = unique(meas_idx(:,:),'rows');

% preallocate the coordinate table
meas_coords = zeros([size(coords,3) 3 size(meas_idx,1)]);

% fill in the coordinates for the chosen beads
for z = 1:size(meas_idx,1)
    for q = 1:size(coords,3)
        meas_coords(q,1:3,z) = coords(meas_idx(z,2)+1,1:3,q);
    end
end

% calculate the mean positions
for z = 1:size(meas_coords,3)
     meas_idx(z,3:5) = mean(meas_coords(:,1:3,z),1);
end

% preallocate the table for the recolored beads
MT_coords = zeros([size(MT_idx,1) 3]);

% log the MT coordinates
for z = 1:size(MT_idx,1)
    MT_coords(z,1:3) = coords(MT_idx(z,1)+1,1:3,1);
end

output.meas_coords = meas_coords;
output.meas_idx = meas_idx;
output.MT_coords = MT_coords;
