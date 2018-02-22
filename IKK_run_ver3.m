function [IKK_measure,FWHM,counts] = IKK_run_ver3(data_cat,reference_spot,red_dist_lim,green_dist_lim)
%% Written by Brandon Friedman
%% User Instructions
% 1) The purpose of this code is to extract measurements on the
% distance between spots (both gaussian and brightest pixel), the positions
% of the red and green spots relative to each other, the FHWM of individual
% spots, the FWHM of an aggretate spot image, and a count of the number of
% cells that fall within chosen parameters.
% 2) The first input is the data_cat table, which is the output for the
% IKK_matfile_combine code.
% 3) The second input is the reference_spot, which denotes the spot about
% which the image should be rotated. This is often the spot that defines
% the spindle access, and this input should be either 'RFP' or 'GFP'
% 4) The parameters are listed below, with descriptions, and mechanisms to
% maintain appropriate values for all parameters are included in the code

%% Parameters
start_pos = 2; % starting row in the data_full file (minimum value is 2)
tilt_lim = 600; % number of nm allowed between chosen z-stacks
gauss_fit_dim = 7; % dimensons of box used for gaussian fitting
gauss_fit_tol = 2; % how many pixels away the new gaussian fitting coordinates can be from the bp
bkd_calc_dim = [7 9]; % dimensions of the box used to calculate the background
gof_thresh = 0.80; % the lower limit for goodness of fit for gaussian fitting
padding = 15; % size of the padding for the logical rotation
mask_hy = 15; % height of the mask for the Y-FWHM measurement
mask_wy = 7; % width of the mask for the Y-FWHM measurement
im_footprint = 0; % binary to turn on footprinting
FWHM.GFP_footprint = {}; % holds the masked GFP footprint for a test
FWHM.R_y = []; % stores the FHWM values
FWHM.G_y = []; % stores the FHWM values
FWHM.R_im_cat = []; % used to make an aggregate FWHM measurement for the RFP
FWHM.G_im_cat = []; % used to make an aggregate FWHM measurement for the GFP

%% Initialization
IKK_measure.meas_bp = {'r1_new','r2_new','g1_new','g2_new','Step Size (nm)','Pixel Size (nm)','GFP_image_stack','RFP_image_stack','gg_dist','rr_dist','GFP1_bkd','GFP2_bkd','RFP1_bkd','RFP2_bkd'};
IKK_measure.meas_gauss = {'r1_new','r2_new','g1_new','g2_new','Step Size (nm)','Pixel Size (nm)','GFP_image_stack','RFP_image_stack','gg_dist','rr_dist'};
counts.bp_success = 0; % count of cells where the bp measurement was successful
counts.bp_fail = 0; % count of cells where the bp measurement failed
counts.gauss_fit_success = 0; % count of cells were gaussian fitting was successful
counts.gauss_fit_fail = 0; % count of cells were gaussian fitting failed

%% Error Messages
if mod(gauss_fit_dim,2) == 0 || mod(mask_hy,2) == 0
    error('gauss_fit_dim and mask_hy must be an odd number');
end

% loop through all the cells
for z = start_pos:size(data_cat.raw_data,1)
    
    % collect the coordinates for all four spots
    r1 = data_cat.raw_data{z,1};
    r2 = data_cat.raw_data{z,2};
    g1 = data_cat.raw_data{z,3};
    g2 = data_cat.raw_data{z,4};
    
    % check for empty cells
    if size(r1,1) ~= 0 && size(r2,1) ~= 0 && size(g1,1) ~= 0 && size(g2,1) ~= 0
        
        % if the spots were misclicked, swap them
        if distance_between_IKK_ver3(r1,g1) <= distance_between_IKK_ver3(r1,g2)
        else
            g1 = data_cat.raw_data{z,4};
            g2 = data_cat.raw_data{z,3};
        end
        
        % obtain pixel size and step size for limits
        pixel_size = data_cat.raw_data{z,6};
        step_size = data_cat.raw_data{z,5};
        plane_separation = tilt_lim/step_size;
        
        % exclude all cells that do not fit the chosen criteria for measurement
        [limits, rr_dist, gg_dist] = IKK_limits_ver3(r1,r2,g1,g2,red_dist_lim,green_dist_lim,plane_separation,pixel_size);
        
        % get cropped data coords
        r1_crop = data_cat.cropped_data{z,1};
        r2_crop = data_cat.cropped_data{z,2};
        g1_crop = data_cat.cropped_data{z,3};
        g2_crop = data_cat.cropped_data{z,4};
        
        % if the spots were misclicked, swap them
        if distance_between_IKK_ver3(r1_crop,g1_crop) <= distance_between_IKK_ver3(r1_crop,g2_crop)
        else
            g1_crop = data_cat.cropped_data{z,4};
            g2_crop = data_cat.cropped_data{z,3};
        end
        
        % pull the planes of interest for the four spots out
        RFP1 = data_cat.cropped_data{z,8}(:,:,r1_crop(3));
        RFP2 = data_cat.cropped_data{z,8}(:,:,r2_crop(3));
        GFP1 = data_cat.cropped_data{z,7}(:,:,g1_crop(3));
        GFP2 = data_cat.cropped_data{z,7}(:,:,g2_crop(3));
        
        % crop the images to an NxN matrix
        RFP1_a = crop_image_IKK_ver3(RFP1, r1_crop, (gauss_fit_dim-1)/2, (gauss_fit_dim-1)/2);
        RFP2_a = crop_image_IKK_ver3(RFP2, r2_crop, (gauss_fit_dim-1)/2, (gauss_fit_dim-1)/2);
        GFP1_a = crop_image_IKK_ver3(GFP1, g1_crop, (gauss_fit_dim-1)/2, (gauss_fit_dim-1)/2);
        GFP2_a = crop_image_IKK_ver3(GFP2, g2_crop, (gauss_fit_dim-1)/2, (gauss_fit_dim-1)/2);
        
        % find the maximum intensity of each image and make sure it's equal to the intensity of the middle pixel
        max_check_r1 = max(RFP1_a(:)) ~= r1(4);
        max_check_r2 = max(RFP2_a(:)) ~= r2(4);
        max_check_g1 = max(GFP1_a(:)) ~= g1(4);
        max_check_g2 = max(GFP2_a(:)) ~= g2(4);
        max_check = max_check_r1 + max_check_r2 + max_check_g1 + max_check_g2;
        
        % if the cell falls within the limits and is the brighest local pixel, then continue
        if limits + max_check == 0
            % make the smaller region for the background measurement
            RFP1_bkd_s = crop_image_IKK_ver3(RFP1, r1_crop, (bkd_calc_dim(1)-1)/2, (bkd_calc_dim(1)-1)/2);
            RFP2_bkd_s = crop_image_IKK_ver3(RFP2, r2_crop, (bkd_calc_dim(1)-1)/2, (bkd_calc_dim(1)-1)/2);
            GFP1_bkd_s = crop_image_IKK_ver3(GFP1, g1_crop, (bkd_calc_dim(1)-1)/2, (bkd_calc_dim(1)-1)/2);
            GFP2_bkd_s = crop_image_IKK_ver3(GFP2, g2_crop, (bkd_calc_dim(1)-1)/2, (bkd_calc_dim(1)-1)/2);
            
            % make the larger region for the background measurement
            RFP1_bkd_l = crop_image_IKK_ver3(RFP1, r1_crop, (bkd_calc_dim(2)-1)/2, (bkd_calc_dim(2)-1)/2);
            RFP2_bkd_l = crop_image_IKK_ver3(RFP2, r2_crop, (bkd_calc_dim(2)-1)/2, (bkd_calc_dim(2)-1)/2);
            GFP1_bkd_l = crop_image_IKK_ver3(GFP1, g1_crop, (bkd_calc_dim(2)-1)/2, (bkd_calc_dim(2)-1)/2);
            GFP2_bkd_l = crop_image_IKK_ver3(GFP2, g2_crop, (bkd_calc_dim(2)-1)/2, (bkd_calc_dim(2)-1)/2);
            
            % calculate the backgrounds
            RFP1_bkd = (sum(RFP1_bkd_l(:))-sum(RFP1_bkd_s(:)))/(bkd_calc_dim(2)^2-bkd_calc_dim(1)^2);
            RFP2_bkd = (sum(RFP2_bkd_l(:))-sum(RFP2_bkd_s(:)))/(bkd_calc_dim(2)^2-bkd_calc_dim(1)^2);
            GFP1_bkd = (sum(GFP1_bkd_l(:))-sum(GFP1_bkd_s(:)))/(bkd_calc_dim(2)^2-bkd_calc_dim(1)^2);
            GFP2_bkd = (sum(GFP2_bkd_l(:))-sum(GFP2_bkd_s(:)))/(bkd_calc_dim(2)^2-bkd_calc_dim(1)^2);
            
            %% RR and GG Logged for Brightest Pixel that Meet Limit Requirements
            meas_bp_size = size(IKK_measure.meas_bp,1);
            IKK_measure.meas_bp(meas_bp_size+1,1:8) = data_cat.cropped_data(z,1:8);
            IKK_measure.meas_bp{meas_bp_size+1,9} = gg_dist;
            IKK_measure.meas_bp{meas_bp_size+1,10} = rr_dist;
            IKK_measure.meas_bp{meas_bp_size+1,11} = GFP1_bkd;
            IKK_measure.meas_bp{meas_bp_size+1,12} = GFP2_bkd;
            IKK_measure.meas_bp{meas_bp_size+1,13} = RFP1_bkd;
            IKK_measure.meas_bp{meas_bp_size+1,14} = RFP2_bkd;
            
            % count and display the successful recording
            counts.bp_success = counts.bp_success + 1;
            disp('BP measurements were successfully recorded');
            
            %% Gaussian Fitting for New Fitted Coordinates
            % use gaussian curve-fitting to calculate the new gaussian-fitted coords
            try
                [r1_g_un] = gauss_center_calc_IKK_ver_3(RFP1_a,gof_thresh);
                [r2_g_un] = gauss_center_calc_IKK_ver_3(RFP2_a,gof_thresh);
                [g1_g_un] = gauss_center_calc_IKK_ver_3(GFP1_a,gof_thresh);
                [g2_g_un] = gauss_center_calc_IKK_ver_3(GFP2_a,gof_thresh);
            catch
                r1_g_un = [1000000 1000000];
                r2_g_un = [1000000 1000000];
                g1_g_un = [1000000 1000000];
                g2_g_un = [1000000 1000000];
            end
            
            % check the distances between the new and old coords
            [g_fit_check_r1] = distance_between_IKK_ver3(r1_g_un,[((gauss_fit_dim+1)/2) ((gauss_fit_dim+1)/2)]);
            [g_fit_check_r2] = distance_between_IKK_ver3(r2_g_un,[((gauss_fit_dim+1)/2) ((gauss_fit_dim+1)/2)]);
            [g_fit_check_g1] = distance_between_IKK_ver3(g1_g_un,[((gauss_fit_dim+1)/2) ((gauss_fit_dim+1)/2)]);
            [g_fit_check_g2] = distance_between_IKK_ver3(g2_g_un,[((gauss_fit_dim+1)/2) ((gauss_fit_dim+1)/2)]);
            
            % check to make sure there were no errors that skewed the new center too far
            if max([g_fit_check_r1,g_fit_check_r2,g_fit_check_g1,g_fit_check_g2]) <= gauss_fit_tol
                
                % convert the coordinates back to make the new cropped image coords
                r1_g = (r1_g_un-((gauss_fit_dim+1)/2)) + r1_crop(1:2);
                r2_g = (r2_g_un-((gauss_fit_dim+1)/2)) + r2_crop(1:2);
                g1_g = (g1_g_un-((gauss_fit_dim+1)/2)) + g1_crop(1:2);
                g2_g = (g2_g_un-((gauss_fit_dim+1)/2)) + g2_crop(1:2);
                
                % add the plane and intensity to the new gaussian coordinates
                r1_g(3:4) = r1_crop(3:4);
                r2_g(3:4) = r2_crop(3:4);
                g1_g(3:4) = g1_crop(3:4);
                g2_g(3:4) = g2_crop(3:4);
                
                % record the information in a new cell array
                meas_gauss_size = size(IKK_measure.meas_gauss,1);
                IKK_measure.meas_gauss{meas_gauss_size+1,1} = r1_g;
                IKK_measure.meas_gauss{meas_gauss_size+1,2} = r2_g;
                IKK_measure.meas_gauss{meas_gauss_size+1,3} = g1_g;
                IKK_measure.meas_gauss{meas_gauss_size+1,4} = g2_g;
                IKK_measure.meas_gauss(meas_gauss_size+1,5:8) = data_cat.cropped_data(z,5:8);
                IKK_measure.meas_gauss{meas_gauss_size+1,9} = distance_between_IKK_ver3(g1_g,g2_g) * pixel_size;
                IKK_measure.meas_gauss{meas_gauss_size+1,10} = distance_between_IKK_ver3(r1_g,r2_g) * pixel_size;
                % count and display the successful recording
                counts.gauss_fit_success = counts.gauss_fit_success + 1;
                disp('Gaussian measurements were successfully recorded');
            else
                % gaussian fitting wasn't accurate enough
                counts.gauss_fit_fail = counts.gauss_fit_fail + 1;
                disp('Could not fit a gaussian to this cell.');
            end
            
            %% Calculate the FWHM in Y
            
            % pad the data in preparation for the logical rotation
            r1_pad = r1_crop(1:2) + padding;
            r2_pad = r2_crop(1:2) + padding;
            g1_pad = g1_crop(1:2) + padding;
            g2_pad = g2_crop(1:2) + padding;
            
            % pad the image files in preparation for the logical rotation
            RFP1_pad = padarray(RFP1, [padding padding]);
            RFP2_pad = padarray(RFP2, [padding padding]);
            GFP1_pad = padarray(GFP1, [padding padding]);
            GFP2_pad = padarray(GFP2, [padding padding]);
            
            % find the angle between the two RFP spots
            if strcmp(reference_spot,'RFP') == 1
                theta = angle_between_IKK_ver3(r1, r2);
            elseif strcmp(reference_spot,'GFP') == 1
                theta = angle_between_IKK_ver3(g1, g2);
            else
                error('reference_spot must either be "GFP" or "RFP"');
            end
            
            % crop the images for the binary mask
            RFP1_b = crop_image_IKK_ver3(RFP1_pad, r1_pad, (mask_hy+(2*padding)-1)/2, (mask_wy+(2*padding)-1)/2);
            RFP2_b = crop_image_IKK_ver3(RFP2_pad, r2_pad, (mask_hy+(2*padding)-1)/2, (mask_wy+(2*padding)-1)/2);
            GFP1_b = crop_image_IKK_ver3(GFP1_pad, g1_pad, (mask_hy+(2*padding)-1)/2, (mask_wy+(2*padding)-1)/2);
            GFP2_b = crop_image_IKK_ver3(GFP2_pad, g2_pad, (mask_hy+(2*padding)-1)/2, (mask_wy+(2*padding)-1)/2);
            
            % make the mask for the spots
            mask = mask_maker_IKK_ver3(mask_wy, mask_hy, padding, -theta);
            
            % create a footprint of the masks on the GFP image
            if im_footprint == 1
                % create white footprints of the images to be added to the GFP image
                g1_fp = mask_maker_footprint_IKK_ver3(mask, g1_pad, size(GFP1_pad,1), size(GFP1_pad,2));
                g2_fp = mask_maker_footprint_IKK_ver3(mask, g2_pad, size(GFP2_pad,1), size(GFP2_pad,2));
                
                % apply the first footprint
                GFP1_test = GFP1_pad .* g1_fp;
                GFP1_test(GFP1_test==0) = max(GFP1_test(:));
                
                % apply the second footprint
                GFP_both_test = GFP1_test .* g2_fp;
                GFP_both_test(GFP_both_test==0) = max(GFP_both_test(:));
                
                % log the footprint in a matrix
                footprint_size = size(FWHM.GFP_footprint,1);
                FWHM.GFP_footprint{footprint_size+1,1} = GFP_both_test;
                
            else
            end
            
            % apply the mask to the images to isolate the spot
            RFP1_c = RFP1_b .* mask;
            RFP2_c = RFP2_b .* mask;
            GFP1_c = GFP1_b .* mask;
            GFP2_c = GFP2_b .* mask;
            
            % rotate the spots to make them horizontal
            RFP1_d = imrotate(RFP1_c,theta,'crop');
            RFP2_d = imrotate(RFP2_c,theta,'crop');
            GFP1_d = imrotate(GFP1_c,theta,'crop');
            GFP2_d = imrotate(GFP2_c,theta,'crop');
            
            clearvars RFP1_new RFP2_new GFP1_new GFP2_new;
            
            % crop the images to remove the padding
            RFP1_new(:,:) = RFP1_d((1+padding):(mask_hy+padding),(1+padding):(mask_wy+padding));
            RFP2_new(:,:) = RFP2_d((1+padding):(mask_hy+padding),(1+padding):(mask_wy+padding));
            GFP1_new(:,:) = GFP1_d((1+padding):(mask_hy+padding),(1+padding):(mask_wy+padding));
            GFP2_new(:,:) = GFP2_d((1+padding):(mask_hy+padding),(1+padding):(mask_wy+padding));
            
            % assign all zeros to the minimum non-zero value
            RFP1_new(RFP1_new==0) = min(RFP1_b(RFP1_b>0));
            RFP2_new(RFP2_new==0) = min(RFP2_b(RFP2_b>0));
            GFP1_new(GFP1_new==0) = min(GFP1_b(GFP1_b>0));
            GFP2_new(GFP2_new==0) = min(GFP2_b(GFP2_b>0));
            
            % calculate the FWHM
            [R1_FWHM_y] = FWHM_calculate_y_IKK_ver3(RFP1_new,size(RFP1_new,1),gof_thresh);
            [R2_FWHM_y] = FWHM_calculate_y_IKK_ver3(RFP2_new,size(RFP2_new,1),gof_thresh);
            [G1_FWHM_y] = FWHM_calculate_y_IKK_ver3(GFP1_new,size(GFP1_new,1),gof_thresh);
            [G2_FWHM_y] = FWHM_calculate_y_IKK_ver3(GFP2_new,size(GFP2_new,1),gof_thresh);
            
            % store the values
            if min([R1_FWHM_y R2_FWHM_y G1_FWHM_y G2_FWHM_y]) > 0
                % get the size of the stored FWHM values
                FWHM_size = size(FWHM.R_y,1);
                % log the values in the table
                FWHM.R_y(FWHM_size+1,1) = R1_FWHM_y;
                FWHM.R_y(FWHM_size+2,1) = R2_FWHM_y;
                FWHM.G_y(FWHM_size+1,1) = G1_FWHM_y;
                FWHM.G_y(FWHM_size+2,1) = G2_FWHM_y;
                
                % concatinate the images for an aggregate analysis
                FWHM.R_im_cat = cat(3,FWHM.R_im_cat,RFP1_new);
                FWHM.R_im_cat = cat(3,FWHM.R_im_cat,fliplr(RFP2_new));
                FWHM.G_im_cat = cat(3,FWHM.G_im_cat,GFP1_new);
                FWHM.G_im_cat = cat(3,FWHM.G_im_cat,fliplr(GFP2_new));
            else
            end
            
        else
            % cell was excluded due to falling outside limits, or pixel-selection error
            counts.bp_fail = counts.bp_fail + 1;
            disp('Cell excluded due to limits.');
        end
    else
    end
end

FWHM.R_mean = mean(FWHM.R_y);
FWHM.G_mean = mean(FWHM.G_y);

% if size(FWHM.R_im_cat,1) > 0
%     % sum the images to get a composite
%     R_agg_im = sum(FWHM.R_im_cat,3);
%     G_agg_im = sum(FWHM.G_im_cat,3);
%
%     % mirror the RFP
%     R_agg_im((size(R_agg_im,1)+1)/2,:) = R_agg_im((size(R_agg_im,1)+1)/2,:)*2;
%     R_agg_im(1:(((size(R_agg_im,1)+1)/2)-1),:) = R_agg_im(1:(((size(R_agg_im,1)+1)/2)-1),:) + flipud(R_agg_im((((size(R_agg_im,1)+1)/2)+1):size(R_agg_im,1),:));
%     R_agg_im((((size(R_agg_im,1)+1)/2)+1):size(R_agg_im,1),:) = flipud(R_agg_im(1:(((size(R_agg_im,1)+1)/2)-1),:));
%
%     % mirror the GFP
%     G_agg_im((size(G_agg_im,1)+1)/2,:) = G_agg_im((size(G_agg_im,1)+1)/2,:)*2;
%     G_agg_im(1:(((size(G_agg_im,1)+1)/2)-1),:) = G_agg_im(1:(((size(G_agg_im,1)+1)/2)-1),:) + flipud(G_agg_im((((size(G_agg_im,1)+1)/2)+1):size(G_agg_im,1),:));
%     G_agg_im((((size(G_agg_im,1)+1)/2)+1):size(G_agg_im,1),:) = flipud(G_agg_im(1:(((size(G_agg_im,1)+1)/2)-1),:));
%
%     % calculate the aggregate FHWM values
%     [FWHM.R_agg_y] = FWHM_calculate_y_IKK_ver3(R_agg_im,size(R_agg_im,1),gof_thresh);
%     [FWHM.G_agg_y] = FWHM_calculate_y_IKK_ver3(G_agg_im,size(G_agg_im,1),gof_thresh);
%
% else
% end