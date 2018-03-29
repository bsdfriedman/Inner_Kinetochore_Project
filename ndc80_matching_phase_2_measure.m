function [ndc80_phase2_output] = ndc80_matching_phase_2_measure(tiff_file_dir)

% parameters
brightest_plane = 4; % this is the plane with the brightest pixel
gauss_fit_dim = 15; % size of the matrix used for gaussian fitting
gauss_fit_dim_h = 15; % height of the matrix used for FHWM gaussian fitting
gauss_fit_dim_w = 7; % width of the matrix used for FWHM gaussian fitting
num_iterations = 250; % the number of times I add noise to the raw image
intens_frac = 3; % fraction of the range of intensity values used to make the noise
gof_thresh = 0.8; % threshold for the gaussian fitting to find the FWHM values
pixel_size = 64; % pixel size in nm
blackout_rad = 3; % the pixel radius used to black out the first spot

%% Identify the folders
% go to the correct directory
cd(tiff_file_dir);
% identify the folder names to get the radii
folder_names = dir('tiff_ndc80_theta_r*');
% preallocate
folder_rad = zeros([size(folder_names,1) 1]);
% pull out all the radii numbers
for z = 1:size(folder_names,1)
    rad_split = strsplit(folder_names(z).name,{'_r','.'});
    folder_rad(z,1) = str2double(rad_split(2));
end
% sort the radii numbers and remove the nan value
folder_rad = sort(folder_rad);
folder_rad = folder_rad(~isnan(folder_rad));

% preallocate
ndc80_phase2_output.images_cropped = zeros([gauss_fit_dim_h gauss_fit_dim_w (size(folder_rad,1)+1)*2]);
ndc80_phase2_output.images_cropped_n = zeros([gauss_fit_dim_h gauss_fit_dim_w num_iterations (size(folder_rad,1)+1)*2]);
ndc80_phase2_output.FWHM_table = zeros([(size(folder_rad,1)+1) num_iterations*2]);
ndc80_phase2_output.spot_idx_bp = cell([size(folder_rad,1) num_iterations]);
ndc80_phase2_output.spot_idx_gauss = cell([size(folder_rad,1) num_iterations]);

% loop through each of the folders
for z = 1:size(folder_rad,1)+1
    
    if z < size(folder_rad,1)+1
        % go to the correct folder
        cd(sprintf('%s%stiff_ndc80_theta_r%d',tiff_file_dir,filesep,folder_rad(z,1)));
    else
        % go to the correct folder
        cd(sprintf('%s%stiff_ndc80_theta_rndc80',tiff_file_dir,filesep));
    end
    
    % identify the tif file and pull out the name
    tiff_files = dir('*.tif');
    tiff_filename = tiff_files.name;
    
    % pull the info from the tif file
    info = imfinfo(tiff_filename);
    
    % pick out the fourth plane
    clearvars im_plane
    im_plane(:,:)= double(imread(tiff_filename, brightest_plane, 'Info', info));
    
    % if this is the first runthrough
    if z == 1
        % preallocate the folder to store noisy full images
        ndc80_phase2_output.images_full = zeros([size(im_plane,1) size(im_plane,1) size(folder_rad,1)+1]);
        ndc80_phase2_output.images_full_n = zeros([size(im_plane,1) size(im_plane,1) num_iterations size(folder_rad,1)+1]);
    else
    end
    
    % find the intensity of the brightest pixel
    pix_max_1 = max(im_plane(:));
    
    % find the coordinates of the brightest pixel
    [pix_max_rows_1, pix_max_cols_1] = find(im_plane == pix_max_1);
    
    % pick the first index of the brightest pixel, in case there are multiple
    spot_idx_1 = [pix_max_cols_1(1),pix_max_rows_1(1)];
    
    % black out the first spot
    im_plane_blackout = im_plane;
    im_plane_blackout(spot_idx_1(2)-blackout_rad:spot_idx_1(2)+blackout_rad,spot_idx_1(1)-blackout_rad:spot_idx_1(1)+blackout_rad) = min(im_plane(:));
    
    % find the intensity of the brightest pixel for the second spot
    pix_max_2 = max(im_plane_blackout(:));
    
    % find the coordinates of the brightest pixel for the second spot
    [pix_max_rows_2, pix_max_cols_2] = find(im_plane_blackout == pix_max_2);
    
    % pick the first index of the brightest pixel for the second spot, in case there are multiple
    spot_idx_2 = [pix_max_cols_2(1),pix_max_rows_2(1)];
    
    % crop the image down to a smaller region to find the FWHM
    im_plane_c1 = crop_image_IKK_ver3(im_plane, spot_idx_1, (gauss_fit_dim_h-1)/2, (gauss_fit_dim_w-1)/2);
    im_plane_c2 = crop_image_IKK_ver3(im_plane, spot_idx_2, (gauss_fit_dim_h-1)/2, (gauss_fit_dim_w-1)/2);
    
    % log the cropped images (without noise)
    ndc80_phase2_output.images_full(:,:,z) = im_plane;
    ndc80_phase2_output.images_cropped(:,:,z*2-1) = im_plane_c1;
    ndc80_phase2_output.images_cropped(:,:,z*2) = im_plane_c2;
    
    % take each image and add noise several times to
    for y = 1:num_iterations
        %% Add noise for the FWHM
        noise_level = round((max(im_plane_c1(:))-min(im_plane_c1(:)))/intens_frac);
        
        % create an array of noise to add to the image
        noise = randi([0 noise_level],gauss_fit_dim_h,gauss_fit_dim_w);
        
        % add the noise to the image
        im_plane_n1 = im_plane_c1 + noise;
        im_plane_n2 = im_plane_c2 + noise;
        
        % save the noisy image
        ndc80_phase2_output.images_cropped_n(:,:,y,z*2-1) = im_plane_n1;
        ndc80_phase2_output.images_cropped_n(:,:,y,z*2) = im_plane_n2;
        
        %% Calculate the FWHM
        % calculate the FWHM in y
        [FWHM_val_1] = FWHM_calculate_y_IKK_ver3(im_plane_n1,gauss_fit_dim_h,gof_thresh);
        [FWHM_val_2] = FWHM_calculate_y_IKK_ver3(im_plane_n2,gauss_fit_dim_h,gof_thresh);
        
        % log the FWHM value
        ndc80_phase2_output.FWHM_table(z,y*2-1) = FWHM_val_1*pixel_size;
        ndc80_phase2_output.FWHM_table(z,y*2) = FWHM_val_2*pixel_size;
        
        %% Add noise for the distance measurements
        %% Add noise for the FWHM
        % create an array of noise to add to the image
        noise_f = randi([0 noise_level],size(im_plane,1),size(im_plane,1));
        
        % add the noise to the image
        im_plane_n_f = im_plane + noise_f;
        
        % save the noisy image
        ndc80_phase2_output.images_full_n(:,:,y,z) = im_plane_n_f;
        
        %% Find the BP coords
        % find the intensity of the new brightest pixel
        pix_max_n1 = max(im_plane_n_f(:));
        
        % find the coordinates of the new brightest pixel
        [pix_max_rows_n1, pix_max_cols_n1] = find(im_plane_n_f == pix_max_n1);
        
        % log this as the spot index
        spot_idx_n1 = [pix_max_cols_n1(1), pix_max_rows_n1(1)];
        
        % black out the first spot
        im_plane_blackout_n = im_plane_n_f;
        im_plane_blackout_n(spot_idx_n1(2)-blackout_rad:spot_idx_n1(2)+blackout_rad,spot_idx_n1(1)-blackout_rad:spot_idx_n1(1)+blackout_rad) = min(im_plane_n_f(:));
        
        % find the intensity of the brightest pixel for the second spot
        pix_max_n2 = max(im_plane_blackout_n(:));
        
        % find the coordinates of the brightest pixel for the second spot
        [pix_max_rows_n2, pix_max_cols_n2] = find(im_plane_blackout_n == pix_max_n2);
        
        % pick the first index of the brightest pixel for the second spot, in case there are multiple
        spot_idx_n2 = [pix_max_cols_n2(1),pix_max_rows_n2(1)];
        
        % pick the first index of the brightest pixel, in case there are multiple
        ndc80_phase2_output.spot_idx_bp{z,y} = [spot_idx_n1,spot_idx_n2];
        
        %% Find the gaussian coords
        % crop the image down to a smaller region to find the FWHM
        im_plane_cn1 = crop_image_IKK_ver3(im_plane_n_f, spot_idx_n1, (gauss_fit_dim-1)/2, (gauss_fit_dim-1)/2);
        im_plane_cn2 = crop_image_IKK_ver3(im_plane_n_f, spot_idx_n2, (gauss_fit_dim-1)/2, (gauss_fit_dim-1)/2);
        
        % pull out the coordinates of the gaussian center (FWHM values not pulled out)
        [idx_gauss_un1] = gauss_center_calc_ndc80_matching(im_plane_cn1,gof_thresh);
        [idx_gauss_un2] = gauss_center_calc_ndc80_matching(im_plane_cn2,gof_thresh);
        
        % convert the coordinates back to make the new cropped image coords
        idx_gauss_n1 = (idx_gauss_un1-((gauss_fit_dim+1)/2)) + spot_idx_n1(1:2);
        idx_gauss_n2 = (idx_gauss_un2-((gauss_fit_dim+1)/2)) + spot_idx_n2(1:2);
        
        % pick the first index of the brightest pixel, in case there are multiple
        ndc80_phase2_output.spot_idx_gauss{z,y} = [idx_gauss_n1,idx_gauss_n2];
    end
end

% create the data table for FWHM
ndc80_phase2_output.FWHM_data = cell([size(folder_rad,1)+2 3]);
ndc80_phase2_output.FWHM_data(1,1:3) = {'Circle Angle','FWHM','n'};
for z = 1:size(folder_rad,1)+1
    if z < size(folder_rad,1)+1
        ndc80_phase2_output.FWHM_data{z+1,1} = folder_rad(z,1);
    else
        ndc80_phase2_output.FWHM_data(z+1,1) = {'ndc80'};
    end
    FWHM_iter = ndc80_phase2_output.FWHM_table(z,:);
    FWHM_iter = FWHM_iter(FWHM_iter>0);
    ndc80_phase2_output.FWHM_data{z+1,2} = mean(FWHM_iter);
    ndc80_phase2_output.FWHM_data{z+1,3} = size(FWHM_iter,2);
end

% create the data table for bp distances
ndc80_phase2_output.distance_data_bp = cell([size(folder_rad,1)+2 7]);
ndc80_phase2_output.distance_data_bp(1,1:7) = {'Circle Angle','ndc80-ame1 dist','ndc80-ame1 xpos','Percent ame1 DNA','Percent ame1 SPB','Percent Same','n'};
for z = 1:size(folder_rad,1)+1
    % add row labels
    if z < size(folder_rad,1)+1
        ndc80_phase2_output.distance_data_bp{z+1,1} = folder_rad(z,1);
    else
        ndc80_phase2_output.distance_data_bp(z+1,1) = {'ndc80'};
    end
    % preallocate
    dist_vals_set = zeros([num_iterations 1]);
    dist_vals_ref = zeros([num_iterations 1]);
    % pull out the coordinates
    for y = 1:num_iterations
        dist_vals_set(y,1) = norm(ndc80_phase2_output.spot_idx_bp{z,y}(1:2)-ndc80_phase2_output.spot_idx_bp{z,y}(3:4))*pixel_size;
        dist_vals_ref(y,1) = norm(ndc80_phase2_output.spot_idx_bp{size(folder_rad,1)+1,y}(1:2)-ndc80_phase2_output.spot_idx_bp{size(folder_rad,1)+1,y}(3:4))*pixel_size;
    end
    % find the spot separation using kk distances
    sep_vals = (dist_vals_ref-dist_vals_set)/2;
    ndc80_phase2_output.distance_data_bp{z+1,2} = mean(abs(sep_vals));
    ndc80_phase2_output.distance_data_bp{z+1,3} = mean(sep_vals);
    ndc80_phase2_output.distance_data_bp{z+1,4} = size(sep_vals(sep_vals>0),1)/size(sep_vals,1);
    ndc80_phase2_output.distance_data_bp{z+1,5} = size(sep_vals(sep_vals<0),1)/size(sep_vals,1);
    ndc80_phase2_output.distance_data_bp{z+1,6} = size(sep_vals(sep_vals==0),1)/size(sep_vals,1);
    ndc80_phase2_output.distance_data_bp{z+1,7} = size(sep_vals,1);
end

% create the data table for gaussian distances
ndc80_phase2_output.distance_data_gauss = cell([size(folder_rad,1)+2 7]);
ndc80_phase2_output.distance_data_gauss(1,1:7) = {'Circle Angle','ndc80-ame1 dist','ndc80-ame1 xpos','Percent ame1 DNA','Percent ame1 SPB','Percent Same','n'};
for z = 1:size(folder_rad,1)+1
    % add row labels
    if z < size(folder_rad,1)+1
        ndc80_phase2_output.distance_data_gauss{z+1,1} = folder_rad(z,1);
    else
        ndc80_phase2_output.distance_data_gauss(z+1,1) = {'ndc80'};
    end
    % preallocate
    dist_vals_set_g = zeros([num_iterations 1]);
    dist_vals_ref_g = zeros([num_iterations 1]);
    gauss_counter = 1; % counter to remove improper fitting
    % pull out the coordinates
    for y = 1:num_iterations
        % make sure the fitting was correct
        if min(ndc80_phase2_output.spot_idx_gauss{z,y})<0 || min(ndc80_phase2_output.spot_idx_gauss{size(folder_rad,1)+1,y})<0
        else
            dist_vals_set_g(gauss_counter,1) = norm(ndc80_phase2_output.spot_idx_gauss{z,y}(1:2)-ndc80_phase2_output.spot_idx_gauss{z,y}(3:4))*pixel_size;
            dist_vals_ref_g(gauss_counter,1) = norm(ndc80_phase2_output.spot_idx_gauss{size(folder_rad,1)+1,y}(1:2)-ndc80_phase2_output.spot_idx_gauss{size(folder_rad,1)+1,y}(3:4))*pixel_size;
            % increase the counter by 1
            gauss_counter = gauss_counter + 1;
        end
    end
    % find the spot separation using kk distances
    sep_vals_g = (dist_vals_ref_g-dist_vals_set_g)/2;
    ndc80_phase2_output.distance_data_gauss{z+1,2} = mean(abs(sep_vals_g));
    ndc80_phase2_output.distance_data_gauss{z+1,3} = mean(sep_vals_g);
    ndc80_phase2_output.distance_data_gauss{z+1,4} = size(sep_vals_g(sep_vals_g>0),1)/size(sep_vals_g,1);
    ndc80_phase2_output.distance_data_gauss{z+1,5} = size(sep_vals_g(sep_vals_g<0),1)/size(sep_vals_g,1);
    ndc80_phase2_output.distance_data_gauss{z+1,6} = size(sep_vals_g(sep_vals_g==0),1)/size(sep_vals_g,1);
    ndc80_phase2_output.distance_data_gauss{z+1,7} = size(sep_vals_g,1);
end
