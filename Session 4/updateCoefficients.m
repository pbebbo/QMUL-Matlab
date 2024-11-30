function coeffs = updateCoefficients(data, params, isRegularized)
    % Update coefficient values using fminsearch
    if isRegularized
        error_func = @(c) regularised_error(data, c, params.lambda);
    else
        error_func = @(c) unregularised_error(data, c);
    end
    
    coeffs = struct();
    coeffs.under = fminsearch(error_func, params.underfit_coeff);
    coeffs.close = fminsearch(error_func, params.close_coeff);
    coeffs.over = fminsearch(error_func, params.overfit_coeff);
end 