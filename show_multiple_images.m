function show_multiple_images(data_cat,num_rows,num_cols)

% initialize
start_pos = 1; % choose the cell in the table that you want to start on

% determine what the ending value for the loop is
if num_rows*num_cols > size(data_cat.cropped_data,1)-1
    loop_max = size(data_cat.cropped_data,1);
elseif num_rows*num_cols == size(data_cat.cropped_data,1)-1
    loop_max = num_rows*num_cols;
end

% loop through the number of
for z = start_pos+1:loop_max
    
end