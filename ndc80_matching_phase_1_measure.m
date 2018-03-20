function [ndc80_match_im_output] = ndc80_matching_phase_1_measure(tiff_file_dir)

% parameters
brightest_plane = 4; % this is the plane with the brightest pixel
gauss_fit_dim = 7; % size of the matrix used for gaussian fitting
num_iterations = 500; % the number of times I add noise to the raw image
intens_frac = 3; % fraction of the range of intensity values used to make the noise
gof_thresh = 0.8; % threshold for the gaussian fitting to find the FWHM values
pixel_size = 64; % pixel size in nm

%% Identify the folders
% go to the correct directory
cd(tiff_file_dir);
% identify the folder names to get the radii
folder_names = dir('tiff_ndc80_circ_r*');
% preallocate
folder_rad = zeros([size(folder_names,1) 1]);
% pull out all the radii numbers
for z = 1:size(folder_names,1)
    rad_split = strsplit(folder_names(z).name,{'_r','.'});
    folder_rad(z,1) = str2double(rad_split(2));
end
% sort the radii numbers
folder_rad = sort(folder_rad);

% preallocate
% ndc80_match_im_output.images_cropped = zeros([gauss_fit_dim gauss_fit_dim size(folder_rad,1)]);
% ndc80_match_im_output.images_cropped_n = zeros([gauss_fit_dim gauss_fit_dim num_iterations size(folder_rad,1)]);
ndc80_match_im_output.FWHM_table = zeros([size(folder_rad,1) num_iterations]);

% loop through each of the folders
for z = 1:size(folder_rad,1)
    % go to the correct folder
    cd(sprintf('%s%stiff_ndc80_circ_r%d',tiff_file_dir,filesep,folder_rad(z,1)));
    
    % identify the tif file and pull out the name
    tiff_files = dir('*.tif');
    tiff_filename = tiff_files.name;
    
    % pull the info from the tif file
    info = imfinfo(tiff_filename);
    
    % pick out the fourth plane
    clearvars im_plane
    im_plane(:,:)= double(imread(tiff_filename, brightest_plane, 'Info', info));
    
    % find the intensity of the brightest pixel
    pix_max_1 = max(im_plane(:));
    
    % find the coordinates of the brightest pixel
    [pix_max_rows, pix_max_cols] = find(im_plane == pix_max_1);
    
    % pick the first index of the brightest pixel, in case there are multiple
    spot_idx_1 = [pix_max_cols(1),pix_max_rows(1)];
    
    % crop the image down to a smaller region to find the FWHM
    im_plane = crop_image_IKK_ver3(im_plane, spot_idx_1, (gauss_fit_dim-1)/2, (gauss_fit_dim-1)/2);
    
    % log the cropped images (without noise)
    ndc80_match_im_output.images_cropped(:,:,z) = im_plane;
    
    % take each image and add noise several times to 
    for y = 1:num_iterations
        
        noise_level = round((max(im_plane(:))-min(im_plane(:)))/intens_frac);
        
        % create an array of noise to add to the image
        noise = randi([0 noise_level],gauss_fit_dim,gauss_fit_dim);
        
        % add the noise to the image
        im_plane_n = im_plane + noise;
        
        % save the noisy image
        ndc80_match_im_output.images_cropped_n(:,:,y,z) = im_plane_n;
        
        % calculate the FWHM in y
        [FWHM_val] = FWHM_calculate_y_IKK_ver3(im_plane_n,gauss_fit_dim,gof_thresh);
        
        % log the FWHM value
        ndc80_match_im_output.FWHM_table(z,y) = FWHM_val*pixel_size;
    end
end

% create the data_table for FWHM
ndc80_match_im_output.FWHM_data = cell([size(folder_rad,1)+1 3]);
ndc80_match_im_output.FWHM_data(1,1:3) = {'Circle Radius','FWHM','n'};
for z = 1:size(folder_rad,1)
    ndc80_match_im_output.FWHM_data{z+1,1} = folder_rad(z,1);
    FWHM_iter = ndc80_match_im_output.FWHM_table(z,:);
    FWHM_iter = FWHM_iter(FWHM_iter>0);
    ndc80_match_im_output.FWHM_data{z+1,2} = mean(FWHM_iter);
    ndc80_match_im_output.FWHM_data{z+1,3} = size(FWHM_iter,2);
end
