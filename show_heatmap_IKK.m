function show_heatmap_IKK(protein_map,pixsize)

% this code displays the heatmap
x = [-7 7]*pixsize; y = x;
imagesc(x,y,protein_map);
colormap hot;