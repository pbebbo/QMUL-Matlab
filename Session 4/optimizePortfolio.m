function [sig, r, w_space] = optimizePortfolio(X, xav, N)
    % Perform portfolio optimization
    % Inputs:
    %   X: normalized returns matrix
    %   xav: mean of each asset's returns
    %   N: number of assets
    % Outputs:
    %   sig: risk values
    %   r: return values
    %   w_space: optimal weights for each return level

    % Calculate covariance matrix
    cov_mat = cov(X);
    e = ones(1,N);
    
    % Define return space
    r = -0.001:0.0001:0.001;
    Nr = length(r);
    sig = zeros(1,Nr);
    w_space = zeros(N,Nr);
    
    % Setup optimization matrices
    A = [cov_mat e' xav'; e 0 0; xav 0 0];
    b = [zeros(N,1); 1; 0];
    
    % Calculate efficient frontier
    for i = 1:Nr
        b(N+2) = r(i);
        w = A\b;
        w0 = w(1:N);
        sig(i) = w0'*cov_mat*w0;
        w_space(:,i) = w0;
    end
end 