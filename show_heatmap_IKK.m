function show_heatmap_IKK(protein_map,pixsize,pos_protein,ref_protein)

% this code displays the heatmap
x = [-7 7]*pixsize; y = x;
heatmap_disp = interp2(protein_map);
imagesc(x,y,heatmap_disp);
colormap hot;

% calcualte the n
cell_n = sum(protein_map(:))/2;

% display the title
title_disp = sprintf('%s vs %s; n=%d',pos_protein,ref_protein,cell_n);
title(title_disp);
xlabel('Position Along Axis (nm)');
ylabel('Position Perpendicular to Axis (nm)');