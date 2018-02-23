function [distance_measure] = mean_stdev_rr_gg(IKK_measure,bp_or_gauss)

% this code calculates the mean and standard deviation of the rr and gg distances

if strcmp(bp_or_gauss,'bp')
    % if the bp coordinates are being used
    distance_measure.rr_mean = mean(cell2mat(IKK_measure.meas_bp(2:end,10)));
    distance_measure.rr_stdev = std(cell2mat(IKK_measure.meas_bp(2:end,10)));
    distance_measure.gg_mean = mean(cell2mat(IKK_measure.meas_bp(2:end,9)));
    distance_measure.gg_stdev = std(cell2mat(IKK_measure.meas_bp(2:end,9)));
elseif strcmp(bp_or_gauss,'gauss')
    % if the gauss coordinates are being used
    distance_measure.rr_mean = mean(cell2mat(IKK_measure.meas_gauss(2:end,10)));
    distance_measure.rr_stdev = std(cell2mat(IKK_measure.meas_gauss(2:end,10)));
    distance_measure.gg_mean = mean(cell2mat(IKK_measure.meas_gauss(2:end,9)));
    distance_measure.gg_stdev = std(cell2mat(IKK_measure.meas_gauss(2:end,9)));
else
    % error if the user types bp_or_gauss in incorrectly
    error('bp_or_gauss needs to be ''bp'' or ''gauss''')
end