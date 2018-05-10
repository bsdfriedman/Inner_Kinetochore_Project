function micro_sim_convert_horiz_cen_n(sim_string)

% make the images from the cen color file
system(sprintf('python.exe ParseBrownian.py -PSF PSF_gain_GFP.txt -out XML_%s_horiz_cen_n -width 100 -height 100 -pixel_size 64 -voxel_depth 300 -focal_planes 7 -noise 10 %s_colors_horiz_cen.txt %s_timepoints.txt',sim_string,sim_string,sim_string));
system(sprintf('python.exe BrownianXMLtoTIFF.py -green -out tiff_%s_horiz_cen_n XML_%s_horiz_cen_n',sim_string,sim_string));

% make the images from the MT color file
system(sprintf('python.exe ParseBrownian.py -PSF PSF_gain_GFP.txt -out XML_%s_%s_n -width 100 -height 100 -pixel_size 64 -voxel_depth 300 -focal_planes 7 -noise 10 %s_colors_%s.txt %s_timepoints.txt',sim_string,'MT',sim_string,'MT',sim_string));
system(sprintf('python.exe BrownianXMLtoTIFF.py -green -out tiff_%s_%s_n XML_%s_%s_n',sim_string,'MT',sim_string,'MT'));
