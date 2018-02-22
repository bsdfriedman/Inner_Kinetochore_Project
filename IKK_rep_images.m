function IKK_rep_images(IKK_measure)

% loop through all images
for z = 2:size(IKK_measure.meas_gauss,1)
    
    % get the coordinates
    r1_crop = IKK_measure.meas_gauss{z,1};
    r2_crop = IKK_measure.meas_gauss{z,2};
    g1_crop = IKK_measure.meas_gauss{z,3};
    g2_crop = IKK_measure.meas_gauss{z,4};
    
    % make sure the spots are on the same plane
    if r1_crop(3) == r2_crop(3) && g1_crop(3) == g2_crop(3)
        % choose the planes
        red_plane = mat2gray(IKK_measure.meas_gauss{z,8}(:,:,r1_crop(3)));
        green_plane = mat2gray(IKK_measure.meas_gauss{z,7}(:,:,g1_crop(3)));
        blue_plane = zeros([size(red_plane,1) size(red_plane,2)]);
        
        % concatinate the planes
        concat_im = cat(3,red_plane,green_plane,blue_plane);
        
        % show the image
        subplot(1,3,1)
        imshow(red_plane);
        image_title = sprintf('RR: %d',IKK_measure.meas_gauss{z,10});
        title(image_title);
        hold on;
        subplot(1,3,2)
        imshow(green_plane);
        image_title = sprintf('GG: %d',IKK_measure.meas_gauss{z,9});
        title(image_title);
        subplot(1,3,3)
        imshow(concat_im);
        hold off;
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        
        waitforbuttonpress;
    else
    end
end