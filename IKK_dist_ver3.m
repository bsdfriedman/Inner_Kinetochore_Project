function [spot_position] = IKK_dist_ver3(IKK_measure,bp_or_gauss,stationary_protein)
% Written by Brandon Friedman

% This code performs calculations on the distance measurements between spots

% assign the distance measurements
if strcmp(bp_or_gauss,'bp') == 1
    rr_dist = cell2mat(IKK_measure.meas_bp(2:end,10));
    gg_dist = cell2mat(IKK_measure.meas_bp(2:end,9));
elseif strcmp(bp_or_gauss,'gauss') == 1
    rr_dist = cell2mat(IKK_measure.meas_gauss(2:end,10));
    gg_dist = cell2mat(IKK_measure.meas_gauss(2:end,9));
else
    error('bp_or_gauss must either be "bp" or "gauss"');
end

% calculate the percent of the time the RFP and GFP is on the outside
spot_position.r_out = sum(rr_dist > gg_dist)/size(rr_dist,1);
spot_position.g_out = sum(gg_dist > rr_dist)/size(rr_dist,1);
spot_position.rg_same = sum(gg_dist == rr_dist)/size(rr_dist,1);

% calculate the kk dist
if strcmp(stationary_protein,'RFP') == 1
    kk_hist = (rr_dist - gg_dist)/2;
elseif strcmp(stationary_protein,'GFP') == 1
    kk_hist = (gg_dist - rr_dist)/2;
else
    error('stationary_protein must either be "RFP" or "GFP"');
end

% calculate distances between cells separately depending on which spot is outside
if strcmp(stationary_protein,'RFP') == 1
    kk_dist_r_out_mean = mean(abs(kk_hist(kk_hist>0)));
    kk_dist_g_out_mean = mean(abs(kk_hist(kk_hist<0)));
    kk_dist_r_out_std = std(abs(kk_hist(kk_hist>0)));
    kk_dist_g_out_std = std(abs(kk_hist(kk_hist<0)));
elseif strcmp(stationary_protein,'GFP') == 1
    kk_dist_r_out_mean = mean(abs(kk_hist(kk_hist<0)));
    kk_dist_g_out_mean = mean(abs(kk_hist(kk_hist>0)));
    kk_dist_r_out_std = std(abs(kk_hist(kk_hist<0)));
    kk_dist_g_out_std = std(abs(kk_hist(kk_hist>0)));
else
    error('stationary_protein must either be "RFP" or "GFP"');
end

% calculate the distances
spot_position.kk_hist = kk_hist;
spot_position.kk_xpos_mean = mean(spot_position.kk_hist);
spot_position.kk_xpos_std = std(spot_position.kk_hist);
spot_position.kk_dist_mean = mean(abs(spot_position.kk_hist));
spot_position.kk_dist_std = std(abs(spot_position.kk_hist));
spot_position.kk_dist_r_out_mean = kk_dist_r_out_mean;
spot_position.kk_dist_r_out_std = kk_dist_r_out_std;
spot_position.kk_dist_g_out_mean = kk_dist_g_out_mean;
spot_position.kk_dist_g_out_std = kk_dist_g_out_std;

