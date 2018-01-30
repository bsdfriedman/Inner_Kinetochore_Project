function IKK_vischeck(data_cat,GFP_or_RFP,show_pixel)
%% Written by Brandon Friedman
%% User Instructions
% 1) This code allows you to look at representative images of cells, by
% clicking through them one at a time
% 2) The code will open the images one at a time, and will move to the next
% image when you click anywhere, until you go through all images or close
% the image window (Clicking the X in the top-right corner of the window);
% the specific order will be spot 1 then spot 2 for each image set, before
% looping to the next image
% 3) It will also show you the specific pixel that was selected using the
% heatmap_GUI
% 4) You should load the the data_full file from any of your condensed
% matlab files, created using the IKK_matfile_combine code
% 5) The data_full variable as an input should remain unchanged
% 6) For the GFP_or_RFP variable, input either 'GFP' if you want to look at
% GFP images or 'RFP' if yoy want to look at RFP images
% 7) If you want to see the pixel that was selected in the image, set
% show_pixel to 1; if you don't need to see this pixel, set show_pixel to 0
% 8) If you want to save an image as a representative image, simply click
% File and Save As and save the file as an image file (e.g. png or jpeg)
% 9) The code will start looping through images from the first row, if you
% want to change this, simply change the start_pos (default value is 2) in the Parameters
% section (note, this value cannot go below 2)

%% Parameters
start_pos = 2; % determines the start point in the data_full file
% this value should be between 2 and the last image in the data set

%% Determination of Images and Coordinates
% looks at the user input to determine whether to show GFP or RFP images
if GFP_or_RFP == 'GFP'
    % if the user input 'GFP', select the GFP image and coordinates
    coords = 3:4;
    image = 7;
elseif GFP_or_RFP == 'RFP'
    % if the user input 'RFP', select the RFP image and coordinates
    coords = 1:2;
    image = 8;
else
    % If the user did not choose either of the above give an error and
    % prompt the user to choose 'GFP' or 
    error('Please choose either GFP or RFP (input the choice as a string)');
end

%% Looping through Images
% loop through all cell rows in the aggregate file
for z = start_pos:size(data_cat.cropped_data,1)
    % loop through the two spots for either GFP or RFP
    for n = coords
        % makes the image fullscreen (although it doesn't work for the
        % first image for some reason
        set(gcf, 'Position', get(0, 'Screensize'));
        % open the figure
        figure(1)
        % loads the image
        imshow(data_cat.cropped_data{z,image}(:,:,data_cat.cropped_data{z,n}(3)),[])
        % parses the x and y coordinates of the chosen pixel
        x = data_cat.cropped_data{z,n}(1);
        y = data_cat.cropped_data{z,n}(2);
        % if the user wants to see the pixel, it plots the chosen pixel on
        % the image
        hold on;
        if show_pixel == 1
            scatter(x,y);
        else
            % if the user doesn't want to see the pixel, it does not plot it
        end
        hold off;
        % waits until the user clicks somewhere to move to the next image
        waitforbuttonpress
    end
end