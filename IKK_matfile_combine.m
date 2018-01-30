function [data_cat] = IKK_matfile_combine(protein_dir,filename_bi)
%% Written by Brandon Friedman
%% User Instructions
% 1) The protein_dir variable is the folder directory in which all images
% and matfiles are held
% 2) If "bi" is not in the image names (i.e. "GFP_001.tif") then filename_bi = 0
% 3) If "bi" is in the image names (i.e. "GFPbi_001.tif") then filename_bi = 1
% 4) If you wish to change the padding size of the representative images,
% change the parameter under the "Initialization" section
% 5) The number of z steps in the image set is set in the "Initialization" section

%% Folder Layout Instructions
% 1) The major folder should have the name of the protein or strain being
% assessed (e.g. "Ame1" or "KBY7043")
% 2) Within this folder, there should be folders of all the individual
% dates given in the format MMDDYYYY (e.g. "01012017", "01022017", etc.) as
% well as an additional folder named "matfiles" that contains all of the
% matlab files from pixel selection
% 3) The folders with the dates (e.g. "01012017") should contain all GFP and RFP image files
% 4) The GFP image should be named "GFP_XXX.tif" or "GFPbi_XXX.tif" (e.g. "GFP_001.tif")
% 5) The RFP image should be named "RFP_XXX.tif" or "RFPbi_XXX.tif" (e.g. "RFP_001.tif")
% 6) Either both the GFP and RFP images should contain bi in the filename,
% or neither of them should
% 7) The "matfiles" folder should contain all matlab files with the naming
% convention "MMDDYYYY_XXX.mat" where MMDDYYYY is the date the images were
% taken, and XXX was the image set used for this matlab file (e.g "01012017_001.mat")

%% Initialization
% set the data_cell as an empty cell array to remove the error
data_cell = {};

% parameters
padding = 15; % padding for cropping
num_z = 7; % number of z stack in image set

% create titles for the output data structure
data_cat.raw_data = {'r1','r2','g1','g2','Step Size (nm)','Pixel Size (nm)','GFP_image_dir','RFP_image_dir'};
data_cat.cropped_data = {'r1_new','r2_new','g1_new','g2_new','Step Size (nm)','Pixel Size (nm)','GFP_image_stack','RFP_image_stack'};

% locate and navigate to the matfiles
matfile_dir = strcat(protein_dir,filesep,'matfiles');
cd(matfile_dir);

% list all of the matlab files in the folder and pull the names out
matfiles = dir('*.mat');
matfiles_list = {matfiles.name};

if filename_bi ~= 0 && filename_bi ~= 1
    error('filename_bi must equal 0 or 1');
else
end

%% Making the Raw Data Portion of the Output File
% loop through the matlab files and begin to aggregate them
for q = 1:size(matfiles,1)
    % open the appropriate file
    load(matfiles_list{1,q});
    
    % check the current size of the structure so new values can be added
    size_old = size(data_cat.raw_data,1);
    
    % load all the coordinates, step size, and pixel size into the new table
    data_cat.raw_data(size_old+1:size_old+size(data_cell,1)-1,1) = data_cell(2:size(data_cell,1),5);
    data_cat.raw_data(size_old+1:size_old+size(data_cell,1)-1,2) = data_cell(2:size(data_cell,1),6);
    data_cat.raw_data(size_old+1:size_old+size(data_cell,1)-1,3) = data_cell(2:size(data_cell,1),1);
    data_cat.raw_data(size_old+1:size_old+size(data_cell,1)-1,4) = data_cell(2:size(data_cell,1),3);
    data_cat.raw_data(size_old+1:size_old+size(data_cell,1)-1,5) = data_cell(2:size(data_cell,1),7);
    data_cat.raw_data(size_old+1:size_old+size(data_cell,1)-1,6) = data_cell(2:size(data_cell,1),8);
    
    % split the string to get the date and number
    matfile_split = regexp(matfiles_list{q},'[_.]','split');
    
    % list the directory to the appropriate image
    if filename_bi == 0 % "bi" is not in the filename
        data_cat.raw_data(size_old+1:size_old+size(data_cell,1)-1,7) = {strcat(protein_dir,filesep,matfile_split{1},filesep,'GFP_',matfile_split{2},'.tif')};
        data_cat.raw_data(size_old+1:size_old+size(data_cell,1)-1,8) = {strcat(protein_dir,filesep,matfile_split{1},filesep,'RFP_',matfile_split{2},'.tif')};
    else % "bi" is in the filename
        data_cat.raw_data(size_old+1:size_old+size(data_cell,1)-1,7) = {strcat(protein_dir,filesep,matfile_split{1},filesep,'GFPbi_',matfile_split{2},'.tif')};
        data_cat.raw_data(size_old+1:size_old+size(data_cell,1)-1,8) = {strcat(protein_dir,filesep,matfile_split{1},filesep,'RFPbi_',matfile_split{2},'.tif')};
    end
end

%% Obtaining Representative Images and New Coordinates
% loop through the matlab files and begin to aggregate them
for q = 1:size(matfiles,1)
    % navigate to the matfiles and open the appropriate file
    cd(matfile_dir);
    load(matfiles_list{1,q});
    
    % check the current size of the structure so new values can be added
    size_old = size(data_cat.cropped_data,1);
    
    % add the step size and pixel data to the output file
    data_cat.cropped_data(size_old+1:size_old+size(data_cell,1)-1,5) = data_cell(2:size(data_cell,1),7);
    data_cat.cropped_data(size_old+1:size_old+size(data_cell,1)-1,6) = data_cell(2:size(data_cell,1),8);
    
    % split the string to get the date and number
    matfile_split = regexp(matfiles_list{q},'[_.]','split');
    
    % loop through all cells to obtain images and new coordinates
    for z = 2:size(data_cell,1)
        % navigate to the appropriate image folder
        cd(strcat(protein_dir,filesep,matfile_split{1}));
        
        try
        % get the coordinates and add the padding to accomodate for padding the image
        r1_pad = data_cell{z,5} + [padding padding 0 0];
        r2_pad = data_cell{z,6} + [padding padding 0 0];
        g1_pad = data_cell{z,1} + [padding padding 0 0];
        g2_pad = data_cell{z,3} + [padding padding 0 0];
        
        % select the region that the rep image will be cropped from
        x_range(1) = min([r1_pad(1) r2_pad(1) g1_pad(1) g2_pad(1)]) - padding;
        x_range(2) = max([r1_pad(1) r2_pad(1) g1_pad(1) g2_pad(1)]) + padding;
        y_range(1) = min([r1_pad(2) r2_pad(2) g1_pad(2) g2_pad(2)]) - padding;
        y_range(2) = max([r1_pad(2) r2_pad(2) g1_pad(2) g2_pad(2)]) + padding;
        
        % reassigns plane value to a number between 1 and 7
        [r1_plane] = plane_conv_IKK_ver3(r1_pad(3));
        [r2_plane] = plane_conv_IKK_ver3(r2_pad(3));
        [g1_plane] = plane_conv_IKK_ver3(g1_pad(3));
        [g2_plane] = plane_conv_IKK_ver3(g2_pad(3));
        
        % get the new coordinates after the crop
        r1_new = [(r1_pad(1)-x_range(1)+1) (r1_pad(2)-y_range(1)+1) r1_plane r1_pad(4)];
        r2_new = [(r2_pad(1)-x_range(1)+1) (r2_pad(2)-y_range(1)+1) r2_plane r2_pad(4)];
        g1_new = [(g1_pad(1)-x_range(1)+1) (g1_pad(2)-y_range(1)+1) g1_plane g1_pad(4)];
        g2_new = [(g2_pad(1)-x_range(1)+1) (g2_pad(2)-y_range(1)+1) g2_plane g2_pad(4)];
        
        % place these new coordinates into the table
        data_cat.cropped_data{size_old+z-1,1} = r1_new;
        data_cat.cropped_data{size_old+z-1,2} = r2_new;
        data_cat.cropped_data{size_old+z-1,3} = g1_new;
        data_cat.cropped_data{size_old+z-1,4} = g2_new;
        
        % choose the appropriate GFP and RFP files
        if filename_bi == 0 % "bi" is not in the filename
            GFP_file = strcat('GFP_',matfile_split{2},'.tif');
            RFP_file = strcat('RFP_',matfile_split{2},'.tif');
        else % "bi" is in the filename
            GFP_file = strcat('GFPbi_',matfile_split{2},'.tif');
            RFP_file = strcat('RFPbi_',matfile_split{2},'.tif');
        end
        
        % get the cropped region to store in the table
        [RFP_region] = choose_region_IKK_ver3(RFP_file,x_range,y_range,r1_pad(3),num_z,padding);
        [GFP_region] = choose_region_IKK_ver3(GFP_file,x_range,y_range,g1_pad(3),num_z,padding);
        
        % log the representative images in the data structure
        data_cat.cropped_data{size_old+z-1,7} = GFP_region;
        data_cat.cropped_data{size_old+z-1,8} = RFP_region;
        
        catch
        end
    end
end