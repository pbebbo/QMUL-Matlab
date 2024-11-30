function params = initializeRegressionParams(isRegularized)
    % Initialize parameters for regression analysis
    params.noise = 0.05;
    params.no_data_pts = 100;
    params.inc = 5;
    params.no_test_pts = 100;
    
    % Initialize coefficients
    params.underfit_coeff = [0,1,1]; % order 2
    params.close_coeff = [0,1,1,1]; % order 3
    params.overfit_coeff = [0,1,1,1,1,1,1,1,1,1]; % order 9
    
    if isRegularized
        params.lambda = exp(-10);
        params.underfit_coeff = params.underfit_coeff';
        params.close_coeff = params.close_coeff';
        params.overfit_coeff = params.overfit_coeff';
    end
end 