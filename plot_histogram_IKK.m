function plot_histogram_IKK(x_or_y,pixel_size,hist_x,hist_y)

if strcmp(x_or_y,'x') == 1
    hist = hist_x;
elseif strcmp(x_or_y,'y') == 1
    hist = hist_y;
else
    error('x_or_y needs to be either "x" or "y"')
end

histogram(hist,'Normalization','probability');
xticklabels([min(hist):max(hist)]*pixel_size);