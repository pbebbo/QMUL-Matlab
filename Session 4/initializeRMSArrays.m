function [x_RMS, y_RMS] = initializeRMSArrays(inc, no_data_pts)
    % Initialize arrays for RMS error tracking
    N_its = length(inc:inc:no_data_pts);
    x_RMS = zeros(1, N_its);
    
    % Structure for RMS values
    y_RMS = struct();
    y_RMS.data_under = zeros(1, N_its);
    y_RMS.data_close = zeros(1, N_its);
    y_RMS.data_over = zeros(1, N_its);
    y_RMS.test_under = zeros(1, N_its);
    y_RMS.test_close = zeros(1, N_its);
    y_RMS.test_over = zeros(1, N_its);
end 