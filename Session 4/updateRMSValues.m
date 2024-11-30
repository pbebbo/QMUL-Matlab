function y_RMS = updateRMSValues(y_RMS, k, data_pts, test_pts, d, coeffs)
    % Update RMS error values for current iteration
    y_RMS.data_under(k) = RMS_error(data_pts(1:d,:), coeffs.under);
    y_RMS.data_close(k) = RMS_error(data_pts(1:d,:), coeffs.close);
    y_RMS.data_over(k) = RMS_error(data_pts(1:d,:), coeffs.over);
    y_RMS.test_under(k) = RMS_error(test_pts, coeffs.under);
    y_RMS.test_close(k) = RMS_error(test_pts, coeffs.close);
    y_RMS.test_over(k) = RMS_error(test_pts, coeffs.over);
end 