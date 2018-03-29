function show_multiple_images(data_cat,num_rows,num_cols,RFP_or_GFP)

% initialize
start_pos = 1; % choose the cell in the table that you want to start on

% determine what the ending value for the loop is
if num_rows*num_cols >= size(data_cat.cropped_data,1)-start_pos
    % if all cells can be included, then throw all of them in there
    loop_max = size(data_cat.cropped_data,1);
else
    % if there are more cells than necessary, then take only the ones you need
    loop_max = (num_rows*num_cols)+1;
end

% loop through the number of
for z = start_pos+1:loop_max
    if strcmp(RFP_or_GFP,'RFP') == 1
        % print the respective RFP image for each line
        subplot(num_rows,num_cols,z-1);
        imshow(data_cat.cropped_data{z,8}(:,:,data_cat.cropped_data{z,1}(3)),[]);
    elseif strcmp(RFP_or_GFP,'GFP') == 1
        % print the respective RFP image for each line
        subplot(num_rows,num_cols,z-1);
        imshow(data_cat.cropped_data{z,7}(:,:,data_cat.cropped_data{z,3}(3)),[]);
    elseif strcmp(RFP_or_GFP,'Both') == 1
        % pull out both images
        RFP_im = mat2gray(data_cat.cropped_data{z,8}(:,:,data_cat.cropped_data{z,1}(3)));
        GFP_im = mat2gray(data_cat.cropped_data{z,7}(:,:,data_cat.cropped_data{z,3}(3)));
        % make a concatinated image
        both_im = cat(3,RFP_im,GFP_im,zeros([size(RFP_im,1) size(RFP_im,2)]));
        % print the image for this row
        subplot(num_rows,num_cols,z-1);
        imshow(both_im);
    else
        error('RFP_or_GFP must be ''RFP'', ''GFP'', or ''Both''')
    end
end