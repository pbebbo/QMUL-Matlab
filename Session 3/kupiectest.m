function [h, pValue, stat] = kupiectest(violations, T, p)
    % Kupiec test for VaR violations
    % violations: number of VaR violations
    % T: total number of observations
    % p: expected violation rate (e.g., 0.05 for 95% VaR)
    
    % Ensure inputs are scalar
    violations = double(violations);
    T = double(T);
    p = double(p);
    
    % Observed violation rate
    p_hat = violations/T;
    
    % Compute likelihood ratio statistic using log properties
    if p_hat == 0
        stat = -2 .* (T .* log(1-p));
    else
        % Use elementwise operations
        stat = -2 .* ((T-violations) .* log(1-p) + violations .* log(p)) + ...
               2 .* ((T-violations) .* log(1-p_hat) + violations .* log(p_hat));
    end
    
    % Handle potential numerical issues
    if any(isnan(stat(:))) || any(isinf(stat(:)))
        warning('Numerical issues in Kupiec test. Setting stat to 0.');
        stat = 0;
    end
    
    % Test against chi-square distribution
    pValue = 1 - chi2cdf(stat, 1);
    h = (pValue < 0.05);  % Reject at 5% significance
end