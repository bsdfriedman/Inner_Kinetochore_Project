function micro_sim_convert(sim_string,beads_inward)

% loop through the script and make the images from the bead color files
for z = 1:size(beads_inward,2)
    sprintf('python.exe ParseBrownian.py -PSF PSF_gain_GFP.txt -out XML_%s_%d -width 100 -height 100 -pixel_size 64 -voxel_depth 300 -focal_planes 7 %s_colors_%d.txt %s_timepoints .txt',sim_string,beads_inward(z),sim_string,beads_inward(z),sim_string)
    system(sprintf('python.exe ParseBrownian.py -PSF PSF_gain_GFP.txt -out XML_%s_%d -width 100 -height 100 -pixel_size 64 -voxel_depth 300 -focal_planes 7 %s_colors_%d.txt %s_timepoints.txt',sim_string,beads_inward(z),sim_string,beads_inward(z),sim_string));
    system(sprintf('python.exe BrownianXMLtoTIFF.py -green -out tiff_%s_%d XML_%s_%d',sim_string,beads_inward(z),sim_string,beads_inward(z)));
end

% make the images from the MT color file
system(sprintf('python.exe ParseBrownian.py -PSF PSF_gain_GFP.txt -out XML_%s_%s -width 100 -height 100 -pixel_size 64 -voxel_depth 300 -focal_planes 7 %s_colors_%s.txt %s_timepoints.txt',sim_string,'MT',sim_string,'MT',sim_string));
system(sprintf('python.exe BrownianXMLtoTIFF.py -green -out tiff_%s_%s XML_%s_%s',sim_string,'MT',sim_string,'MT'));
