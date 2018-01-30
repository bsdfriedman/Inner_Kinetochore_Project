function [hist_x, hist_y] = make_histograms_IKK(protein_map)

% preallocate the histogram terms
hist_x = zeros([sum(protein_map(:)) 1]);
hist_y = zeros([sum(protein_map(:)) 1]);

% sum the protein map in x and y
x_sum = sum(protein_map,1);
y_sum = sum(protein_map,2);

% preset the counts
x_count = 1;
y_count = 1;

% loop through and make the x-dist histogram
for z = 1:size(protein_map,2)
    if x_sum(z) > 0
        for q = 1:x_sum(z)
            hist_x(x_count,1) = z-((size(protein_map,2)+1)/2);
            x_count = x_count + 1;
        end
    else
    end
end

% loop through and make the y-dist histogram
for z = 1:size(protein_map,1)
    if y_sum(z) > 0
        for q = 1:y_sum(z)
            hist_y(y_count,1) = z-((size(protein_map,1)+1)/2);
            y_count = y_count + 1;
        end
    else
    end
end

