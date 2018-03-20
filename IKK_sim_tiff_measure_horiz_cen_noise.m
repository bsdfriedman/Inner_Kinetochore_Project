function [sim_output] = IKK_sim_tiff_measure_horiz_cen_noise(tiff_dir)

% parameters
blackout_rad = 3; % the pixel radius used to black out the first spot
brightest_plane = 4; % the brightest plane in the simulation image
pixel_size = 64; % pixel size in nm
gauss_fit_dim = 9; % size of the matrix used for gaussian fitting
gof_thresh = 0.5; % threshold for goodness of fit, this should be deliberately low

% split the tiff directory to get the tiff string
tiff_dir_split = strsplit(tiff_dir,sprintf('%s',filesep));
tiff_str_split = strsplit(tiff_dir_split{size(tiff_dir_split,2)},'_');
tiff_str = sprintf('%s_%s_%s_%s',tiff_str_split{1},tiff_str_split{2},tiff_str_split{3},tiff_str_split{4});

%% Identify the bp and gauss coordinates, and the KK distances
% loop through all the folders and collect the tiffs you want
for z = 1:2
    if z == 1
        % go to the tiff folder
        cd(sprintf('%s%stiff_IKK_pericentromere_%s_horiz_cen_n',tiff_dir,filesep,tiff_str));
    else
        % go to the tiff folder
        cd(sprintf('%s%stiff_IKK_pericentromere_%s_MT_n',tiff_dir,filesep,tiff_str));
    end
    
    % identify all the tif files
    tiff_files = dir('*.tif');
    if z == 1
        % preallocate
        sim_output.coords_bp = cell([3 7 size(tiff_files,1)]);
        sim_output.coords_gauss = cell([3 4 size(tiff_files,1)]);
        % add labels to the output table
        for h = 1:size(tiff_files,1)
            sim_output.coords_bp(1,:,h) = {'Bead Idx','Spot 1','Spot 2','Image','KK','FWHM 1 [X,Y]','FWHM 2 [X,Y]'};
            sim_output.coords_gauss(1,:,h) = {'Bead Idx','Spot 1','Spot 2','KK'};
            sim_output.coords_bp{2,1,h} = 'Cen';
            sim_output.coords_gauss{2,1,h} = 'Cen';
            sim_output.coords_bp{3,1,h} = 'MT';
            sim_output.coords_gauss{3,1,h} = 'MT';
        end
    else
    end
    %% Identify the BP coords and KK_distances
    % loop through all the tiff files
    for h = 1:size(tiff_files,1)
        % choose the filename and pull the image info out
        tiff_filename = sprintf('IKK_pericentromere_%s_timepoints_%d_0000_G_%d.tif',tiff_str,h,h);
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
        
        % black out the first spot
        im_plane_blackout = im_plane;
        im_plane_blackout(spot_idx_1(2)-blackout_rad:spot_idx_1(2)+blackout_rad,spot_idx_1(1)-blackout_rad:spot_idx_1(1)+blackout_rad) = min(im_plane(:));
        
        % find the intensity of the brightest pixel for the second spot
        pix_max_2 = max(im_plane_blackout(:));
        
        % find the coordinates of the brightest pixel for the second spot
        [pix_max_rows, pix_max_cols] = find(im_plane_blackout == pix_max_2);
        
        % pick the first index of the brightest pixel for the second spot, in case there are multiple
        spot_idx_2 = [pix_max_cols(1),pix_max_rows(1)];
        
        % calculate the KK distance
        kk_dist_bp = norm(spot_idx_1-spot_idx_2)*pixel_size;
        
        % log the coordinates and the intensity into the matrix
        sim_output.coords_bp{z+1,2,h} = [spot_idx_1(1),spot_idx_1(2),brightest_plane,pix_max_1];
        sim_output.coords_bp{z+1,3,h} = [spot_idx_2(1),spot_idx_2(2),brightest_plane,pix_max_2];
        sim_output.coords_bp{z+1,4,h} = im_plane;
        sim_output.coords_bp{z+1,5,h} = kk_dist_bp;
        
        %% Identify the Gauss coords and KK_distances
        % crop the image down to a smaller region to find the gaussian-fitted center
        im_plane_a1 = crop_image_IKK_ver3(im_plane, spot_idx_1, (gauss_fit_dim-1)/2, (gauss_fit_dim-1)/2);
        im_plane_a2 = crop_image_IKK_ver3(im_plane, spot_idx_2, (gauss_fit_dim-1)/2, (gauss_fit_dim-1)/2);
        
        % perform a gaussian fitting and pull out the unfixed coordinates
        [idx_gauss_1_un,FHWM_spot_1] = gauss_sim_fit_IKK_ver_3(im_plane_a1,gof_thresh);
        [idx_gauss_2_un,FHWM_spot_2] = gauss_sim_fit_IKK_ver_3(im_plane_a2,gof_thresh);
        
        % convert the coordinates back to make the new cropped image coords
        idx_gauss_1 = (idx_gauss_1_un-((gauss_fit_dim+1)/2)) + spot_idx_1(1:2);
        idx_gauss_2 = (idx_gauss_2_un-((gauss_fit_dim+1)/2)) + spot_idx_2(1:2);
        
        % calculate the KK_distance
        kk_dist_gauss = norm(idx_gauss_1-idx_gauss_2)*pixel_size;
        
        % log the coordinates and the intensity into the matrix
        sim_output.coords_gauss{z+1,2,h} = [idx_gauss_1(1),idx_gauss_1(2),brightest_plane,pix_max_1];
        sim_output.coords_gauss{z+1,3,h} = [idx_gauss_2(1),idx_gauss_2(2),brightest_plane,pix_max_2];
        sim_output.coords_gauss{z+1,4,h} = kk_dist_gauss;
        sim_output.coords_bp{z+1,6,h} = FHWM_spot_1;
        sim_output.coords_bp{z+1,7,h} = FHWM_spot_2;
    end
end
%% Perform measurements based on the extracted numbers
% preallocate
kk_vals_bp = zeros([size(sim_output.coords_bp,3) 2]);
kk_vals_gauss = zeros([size(sim_output.coords_bp,3) 2]);

% pull out all the kk values
for z = 1:2
    for h = 1:size(sim_output.coords_bp,3)
        kk_vals_bp(h,z) = sim_output.coords_bp{z+1,5,h};
        kk_vals_gauss(h,z) = sim_output.coords_gauss{z+1,4,h};
    end
end

% preallocate
sim_output.kk_values = cell([3 5]);
sim_output.kk_values(1,:) = {'Bead Idx','KK Mean (bp)','KK StDev (bp)','KK Mean (gauss)','KK StDev (gauss)'};
sim_output.kk_values{2,1} = 'Cen';
sim_output.kk_values{3,1} = 'MT';

% add the mean and stdev of the kk values to the table
for z = 1:2
    sim_output.kk_values{z+1,2} = mean(kk_vals_bp(:,z));
    sim_output.kk_values{z+1,3} = std(kk_vals_bp(:,z));
    sim_output.kk_values{z+1,4} = mean(kk_vals_gauss(:,z));
    sim_output.kk_values{z+1,5} = std(kk_vals_gauss(:,z));
end

% preallocate
sim_output.DNA_cen_pos_bp = cell([3 8]);
sim_output.DNA_cen_pos_bp(1,:) = {'Bead Idx','Cen-Dist','Cen-Dist StDev','Cen-Xpos Mean','Cen-Xpos StDev','Percent SPB','Percent DNA','Percent Same'};
sim_output.DNA_cen_pos_bp{2,1} = 'Cen';
sim_output.DNA_cen_pos_bp{3,1} = 'MT';

% find the bp measurements for MT-cen measurements
for z = 1:2
    kk_meas = (kk_vals_bp(:,1)-kk_vals_bp(:,z))/2;
    sim_output.DNA_cen_pos_bp{z+1,2} = mean(abs(kk_meas));
    sim_output.DNA_cen_pos_bp{z+1,3} = std(abs(kk_meas));
    sim_output.DNA_cen_pos_bp{z+1,4} = mean(kk_meas);
    sim_output.DNA_cen_pos_bp{z+1,5} = std(kk_meas);
    sim_output.DNA_cen_pos_bp{z+1,6} = size(kk_meas(kk_meas<0),1)/size(kk_meas,1);
    sim_output.DNA_cen_pos_bp{z+1,7} = size(kk_meas(kk_meas>0),1)/size(kk_meas,1);
    sim_output.DNA_cen_pos_bp{z+1,8} = size(kk_meas(kk_meas==0),1)/size(kk_meas,1);
end

% preallocate
sim_output.DNA_cen_pos_gauss = cell([3 8]);
sim_output.DNA_cen_pos_gauss(1,:) = {'Bead Idx','Cen-Dist','Cen-Dist StDev','Cen-Xpos Mean','Cen-Xpos StDev','Percent SPB','Percent DNA','Percent Same'};
sim_output.DNA_cen_pos_gauss{2,1} = 'Cen';
sim_output.DNA_cen_pos_gauss{3,1} = 'MT';

% find the gauss measurements for DNA-cen measurements
for z = 1:2
    kk_meas = (kk_vals_gauss(:,1)-kk_vals_gauss(:,z))/2;
    sim_output.DNA_cen_pos_gauss{z+1,2} = mean(abs(kk_meas));
    sim_output.DNA_cen_pos_gauss{z+1,3} = std(abs(kk_meas));
    sim_output.DNA_cen_pos_gauss{z+1,4} = mean(kk_meas);
    sim_output.DNA_cen_pos_gauss{z+1,5} = std(kk_meas);
    sim_output.DNA_cen_pos_gauss{z+1,6} = size(kk_meas(kk_meas<0),1)/size(kk_meas,1);
    sim_output.DNA_cen_pos_gauss{z+1,7} = size(kk_meas(kk_meas>0),1)/size(kk_meas,1);
    sim_output.DNA_cen_pos_gauss{z+1,8} = size(kk_meas(kk_meas==0),1)/size(kk_meas,1);
end

% preallocate
sim_output.FHWM = cell([3 5]);
sim_output.FHWM(1,:) = {'Bead Idx','FWHM-X Mean','FHWM-X StDev','FWHM-Y Mean','FHWM-Y StDev'};
sim_output.FHWM{2,1} = 'Cen';
sim_output.FHWM{3,1} = 'MT';

% preallocate
x_FHWM = zeros([size(sim_output.coords_bp,3)*2 2]);
y_FHWM = zeros([size(sim_output.coords_bp,3)*2 2]);

% pull out the FHWM data
for z = 1:2
    for h = 1:size(sim_output.coords_bp,3)
        x_FHWM(h*2-1,z) = sim_output.coords_bp{z+1,6,h}(1);
        x_FHWM(h*2,z) = sim_output.coords_bp{z+1,7,h}(1);
        y_FHWM(h*2-1,z) = sim_output.coords_bp{z+1,6,h}(2);
        y_FHWM(h*2,z) = sim_output.coords_bp{z+1,7,h}(2);
    end
end

% calcualte the means and stdev for FWHM
for z = 1:2
    sim_output.FHWM{z+1,2} = mean(x_FHWM(:,z))*pixel_size;
    sim_output.FHWM{z+1,3} = std(x_FHWM(:,z))*pixel_size;
    sim_output.FHWM{z+1,4} = mean(y_FHWM(:,z))*pixel_size;
    sim_output.FHWM{z+1,5} = std(y_FHWM(:,z))*pixel_size;
end

