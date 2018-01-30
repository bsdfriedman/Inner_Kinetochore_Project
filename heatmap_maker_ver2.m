function [protein_map] = heatmap_maker_ver2(data_cat,stationary_color)

start_pos = 2; % starting row in the data_full file (minimum value is 2)
tilt_lim = 600; % number of nm allowed between chosen z-stacks
red_dist_lim = [0 2000]; % range of acceptable red-red distances
green_dist_lim = [0 2000]; % range of acceptable green-green distances
gauss_fit_dim = 7; % dimensons of box used for gaussian fitting
heatmap_size = 15; % size of the heatmap image
sf = 1; % scaling factor for pixels
plane = zeros([heatmap_size*sf heatmap_size*sf (size(data_cat.raw_data,1)-1)]); % preallocate the planes

if mod(heatmap_size,2) == 0 || mod(sf,2) == 0
    error('Both heatmap_size and sf need to be odd integers');
else
end

% loop through all the cells
for z = start_pos:size(data_cat.raw_data,1)
    
    % collect the coordinates for all four spots
    r1 = data_cat.raw_data{z,1};
    r2 = data_cat.raw_data{z,2};
    g1_t = data_cat.raw_data{z,3};
    g2_t = data_cat.raw_data{z,4};
    
    % make sure the spots are on the correct side
    if distance_between_IKK_ver3(r1,g1_t) < distance_between_IKK_ver3(r1,g2_t)
        g1 = g1_t;
        g2 = g2_t;
    else
        g1 = g2_t;
        g2 = g1_t;
    end
    
    % obtain pixelsize and step size for limits
    pixel_size = data_cat.raw_data{z,6};
    step_size = data_cat.raw_data{z,5};
    plane_separation = tilt_lim/step_size;
    
    % exclude all cells that do not fit the chosen criteria for measurement
    [limits] = IKK_limits_ver3(r1,r2,g1,g2,red_dist_lim,green_dist_lim,plane_separation,pixel_size);
    
    % get cropped data coords
    r1_crop = data_cat.cropped_data{z,1};
    r2_crop = data_cat.cropped_data{z,2};
    g1_crop = data_cat.cropped_data{z,3};
    g2_crop = data_cat.cropped_data{z,4};
    
    % pull the planes of interest for the four spots out
    RFP1 = data_cat.cropped_data{z,8}(:,:,r1_crop(3));
    RFP2 = data_cat.cropped_data{z,8}(:,:,r2_crop(3));
    GFP1 = data_cat.cropped_data{z,7}(:,:,g1_crop(3));
    GFP2 = data_cat.cropped_data{z,7}(:,:,g2_crop(3));
    
    % crop the images to an NxN matrix
    RFP1_a = crop_image_IKK_ver3(RFP1, r1_crop, gauss_fit_dim, gauss_fit_dim);
    RFP2_a = crop_image_IKK_ver3(RFP2, r2_crop, gauss_fit_dim, gauss_fit_dim);
    GFP1_a = crop_image_IKK_ver3(GFP1, g1_crop, gauss_fit_dim, gauss_fit_dim);
    GFP2_a = crop_image_IKK_ver3(GFP2, g2_crop, gauss_fit_dim, gauss_fit_dim);
    
    % find the maximum intensity of each image and make sure it's equal to the intensity of the middle pixel
    max_check_r1 = max(RFP1_a(:)) ~= r1(4);
    max_check_r2 = max(RFP2_a(:)) ~= r2(4);
    max_check_g1 = max(GFP1_a(:)) ~= g1(4);
    max_check_g2 = max(GFP2_a(:)) ~= g2(4);
    max_check = max_check_r1 + max_check_r2 + max_check_g1 + max_check_g2;
    
    % determine which spots are stationary, and which are dynamic, based on the user input
    if strcmp(stationary_color,'GFP')
        s1 = g1(1:2);
        s2 = g2(1:2);
        d1 = r1(1:2);
        d2 = r2(1:2);
    elseif strcmp(stationary_color,'RFP')
        s1 = r1(1:2);
        s2 = r2(1:2);
        d1 = g1(1:2);
        d2 = g2(1:2);
    else
        error('Stationary Color needs to be "GFP" or "RFP"')
    end
    
    % produce the new normalized coordinates
    s2_new = s2-s1;
    d1_new = d1-s1;
    d2_new = d2-s1;
    
    % if the cell falls within the limits and is the brighest local pixel, then continue
    if limits + max_check == 0
        % rotate and round the coordinates
        [theta] = angle_between_IKK_ver3(s1, s2);
        rot_mat = [cos(-theta),-sin(-theta);sin(-theta),cos(-theta)];
        s2_rot = round(sf*(rot_mat*s2_new'))';
        d1_rot = round(sf*(rot_mat*d1_new'))';
        d2_rot = round(sf*(rot_mat*d2_new'))';
        
        % lip the second spot
        d2_fix = d2_rot-s2_rot;
        d2_fix(1) = d2_fix(1)*-1;
        
        % plot these two spots on the heatmap plane
        mid = ((heatmap_size*sf)+1)/2;
        d1_plot = d1_rot + mid;
        d2_plot = d2_fix + mid;
        
        % add the spots to the heatmap
        plane(d1_plot(2),d1_plot(1),z-1) = plane(d1_plot(2),d1_plot(1),z-1) + 1;
        plane(d2_plot(2),d2_plot(1),z-1) = plane(d2_plot(2),d2_plot(1),z-1) + 1;
        
    else
    end
end

% sum the matrices to get the final heatmap
protein_map = sum(plane,3);

% mirror the heatmap
prot_map_temp = zeros([size(protein_map,1) size(protein_map,2)]);
prot_map_temp(((size(protein_map,1)+1)/2),:) = protein_map(((size(protein_map,1)+1)/2),:)*2;
for z = 1:((size(protein_map,1)-1)/2)
    prot_map_temp(z,:) = protein_map(z,:) + protein_map((size(protein_map,1)+1-z),:);
    prot_map_temp((size(protein_map,1)+1-z),:) = protein_map(z,:) + protein_map((size(protein_map,1)+1-z),:);
end

% reassign the protein map
protein_map = prot_map_temp;
