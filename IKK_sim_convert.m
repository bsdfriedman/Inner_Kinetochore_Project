function [output] = IKK_sim_convert(infile,beads_inward,every_n_timesteps)

% The beads_inward input denotes the number of beads inward from the
% centromere you want to be labeled. For example, if it is 1, the two beads
% adjacent to each centromere will be labeled white. You can create
% multiple color files by making a larger array (E.g. [2 5] will make two
% color files with bead 2 labeled in one, and bead 5 labeled in the other).

%% Initialize
% log the time steps that will be used in the calculation
[~,mass_num] = system(sprintf('grep "mass " %s | wc -l',infile));
mass_num = str2double(mass_num);

% log the time steps that will be used in the calculation
[~,times] = system(sprintf('grep "Time " %s | gawk ''!(NR%%%d)''',infile,every_n_timesteps));
time_name = strsplit(times,'\n')';
time_name = time_name(~cellfun('isempty',time_name));

% parameters
mass_mass = 3.38889e-020; % mass of a standard bead
mass_multiplier = 10^10; % pinning mass
ignore = 0; % how many time points you want to ignore

% preallocate
cen_idx = zeros([16 1]);
MT_idx = zeros([16 1]);
coords = zeros([mass_num 3 size(time_name,1)]);

% counters
cen_counter = 1;
time_counter = 1;
MT_counter = 1;

%% Create the timepoint files, log colors, and log coordinates
% open up the file
fid_in = fopen(infile);

% assign the name of the new file that is to be created
file_root_name = strsplit(infile,'.');
filename = sprintf('%s_timepoints.txt',file_root_name{1});
fid_out = fopen(filename,'w');

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
    elseif size(strfind(tline,'Time '),1) ~= 0
        % reprint the line
        fprintf(fid_out,'%s\r\n',tline);
        tline_rec = tline;
        % loop to the next line
        tline = fgetl(fid_in);
        for z = 1:mass_num
            b = str2double(strsplit(tline))*10^6;
            % reprint the line with the new measurements
            fprintf(fid_out,'%0.6f %0.6f %0.6f\r\n',b(1),b(2),b(3));
            if max(strcmp(tline_rec,time_name)) == 1
                % log the coordinates
                coords(z,1:3,time_counter) = str2double(strsplit(tline));
            else
            end
            % loop to the next line
            tline = fgetl(fid_in);
        end
        if max(strcmp(tline_rec,time_name)) == 1
            % increase the time counter by 1
            time_counter = time_counter + 1;
        else
        end
    else
        % loop to the next line
        tline = fgetl(fid_in);
    end
end

% close all the files
fclose('all');

% flip the coordinates so that they are right side up
coords = flipud(coords);

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

% if the cen beads are chosen, remove duplicate values
recolor_idx = unique(recolor_idx,'rows');

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
    
    for z = 1:size(mass_color_list_new,1)
        % print out the color names into the color file
        fprintf(fid_out,'%d\r\n',mass_color_list_new(z,1,q));
    end
    
    % close all the files
    fclose('all');
end

%% THIS CODE IS USED TO TEST THE PROGRAM AND HAS BEEN COMMENTED OUT
% for q = 1:size(beads_inward,2)
%     % open up the file
%     fid_in = fopen(infile);
%
%     % create the colors file
%     filename = sprintf('%s_reprint_%d.out',file_root_name{1},beads_inward(q));
%     fid_out = fopen(filename,'w');
%
%     % assign tline so that the lines can be looped through
%     tline = fgetl(fid_in);
%
%     while ischar(tline) == 1
%         if size(strfind(tline,'MassColors'),1) ~= 0
%             % reprint the line
%             fprintf(fid_out,'%s\r\n',tline);
%             % loop to the next line
%             tline = fgetl(fid_in);
%             for z = 1:size(mass_color_list,1)
%                 % print out the color names into the color file
%                 fprintf(fid_out,'%d\r\n',mass_color_list_new(z,1,q));
%                 % loop to the next line
%                 tline = fgetl(fid_in);
%             end
%         else
%             % reprint the line
%             fprintf(fid_out,'%s\r\n',tline);
%             % loop to the next line
%             tline = fgetl(fid_in);
%         end
%     end
%
%     % close all the files
%     fclose('all');
% end

%% Performing the simulation measurements
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

% create the outputs
output.cen_bead_dist_table = cen_bead_dist_table;
output.cen_bead_xpos_table = cen_bead_xpos_table;
output.MT_cen_xpos_table = MT_cen_xpos_table;
output.cen_bead_dist = cen_bead_dist;
output.cen_bead_xpos = cen_bead_xpos;
output.xpos_percent = xpos_percent;
output.MTpos_percent = MTpos_percent;
