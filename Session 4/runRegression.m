function runRegression(isRegularized)
    % Main function to run regression analysis
    params = initializeRegressionParams(isRegularized);
    
    % Generate data
    data_pts = generate_data(@f, params.noise, params.no_data_pts);
    test_pts = generate_data(@f, params.noise, params.no_test_pts);
    x = linspace(0,1,1000);
    y = f(x);
    
    % Initialize RMS arrays
    [x_RMS, y_RMS] = initializeRMSArrays(params.inc, params.no_data_pts);
    
    % Main loop
    for d = params.inc:params.inc:params.no_data_pts
        % Update coefficients
        coeffs = updateCoefficients(data_pts(1:d,:), params, isRegularized);
        
        % Generate predictions
        y_pred = struct();
        y_pred.under = predictPolynomial(x, coeffs.under);
        y_pred.close = predictPolynomial(x, coeffs.close);
        y_pred.over = predictPolynomial(x, coeffs.over);
        
        % Update RMS values
        k = d/params.inc;
        x_RMS(k) = d;
        y_RMS = updateRMSValues(y_RMS, k, data_pts, test_pts, d, coeffs);
        
        % Plot results
        plotRegressionResults(x, y, data_pts, d, y_pred, x_RMS, y_RMS);
        input('Press Enter to continue...\n');
    end
end 