function IKK_pericen_recolor(infile,beads_inward)

% The beads_inward input denotes the number of beads inward from the
% centromere you want to be labeled. For example, if it is 1, the two beads
% adjacent to each centromere will be labeled white. This can be a single
% bead (e.g. "5" for position 5) or multiple beads (e.g. [2 5] for both
% positions 5 and 2.

% parameters
mass_mass = 3.38889e-020; % mass of a standard bead
mass_multiplier = 10^10; % pinning mass

% preallocate
cen_idx = zeros([16 1]);

% counters
cen_counter = 1;

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
        else
        end
        mass_last = str2double(b{3}); % assigns the last mass value in the file
        % loop to the next line
        tline = fgetl(fid_in);
    elseif size(strfind(tline,'MassColors'),1) ~= 0
        mass_color_list = zeros([mass_last+1 1]);
        if beads_inward > cen_idx(1,1)
            % if the chosen position exceeds the maximum possible value
            error('This position is too far from the centromere')
        else
        end
        for z = 1:mass_last+1
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
        % loop to the next line
        tline = fgetl(fid_in);
        for z = 1:mass_last+1
            b = str2double(strsplit(tline))*10^6;
            % reprint the line with the new measurements
            fprintf(fid_out,'%0.6f %0.6f %0.6f\r\n',b(1),b(2),b(3));
            % loop to the next line
            tline = fgetl(fid_in);
        end
    else
        % loop to the next line
        tline = fgetl(fid_in);
    end
end

% close all the files
fclose('all');

% preallocate the table for the recolored beads
recolor_idx = zeros([32 1 ]);

% choose the beads that are to be recolored
for q = 1:size(beads_inward,2)
    for z = 1:size(cen_idx,1)
        recolor_idx(z*2-1,q) = cen_idx(z,1)-beads_inward(q);
        recolor_idx(z*2,q) = cen_idx(z,1)+beads_inward(q);
    end
end

% if the cen beads are chose, remove duplicate values
recolor_idx = unique(recolor_idx,'rows');

% change the appropriate colors in the mass color list
for q = 1:size(beads_inward,2)
    for z = 1:size(recolor_idx,1)
        mass_color_list(mass_last-recolor_idx(z,q)+1,1) = 4;
    end
end

% create the colors file
filename = sprintf('%s_colors.txt',file_root_name{1});
fid_out = fopen(filename,'w');

for z = 1:size(mass_color_list,1)
    % print out the color names into the color file
    fprintf(fid_out,'%d\r\n',mass_color_list(z,1));
end

% close all the files
fclose('all');

% THE REMAINING CODE IS USED TO TEST THE PROGRAM AND HAS BEEN COMMENTED OUT

% % open up the file
% fid_in = fopen(infile);
% 
% % create the colors file
% filename = sprintf('%s_reprint.out',file_root_name{1});
% fid_out = fopen(filename,'w');
% 
% % assign tline so that the lines can be looped through
% tline = fgetl(fid_in);
% 
% while ischar(tline) == 1
%     if size(strfind(tline,'MassColors'),1) ~= 0
%         % reprint the line
%         fprintf(fid_out,'%s\r\n',tline);
%         % loop to the next line
%         tline = fgetl(fid_in);
%         for z = 1:size(mass_color_list,1)
%             % print out the color names into the color file
%             fprintf(fid_out,'%d\r\n',mass_color_list(z,1));
%             % loop to the next line
%             tline = fgetl(fid_in);
%         end
%     else
%         % reprint the line
%         fprintf(fid_out,'%s\r\n',tline);
%         % loop to the next line
%         tline = fgetl(fid_in);
%     end
% end
% 
% % close all the files
% fclose('all');
