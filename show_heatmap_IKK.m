function show_heatmap_IKK(protein_map,pixsize)

% this code displays the heatmap
x = [-7 7]*pixsize; y = x;
heatmap_disp = interp2(protein_map);
imagesc(x,y,heatmap_disp);
colormap hot;