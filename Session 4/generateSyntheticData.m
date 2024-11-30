function [X, xav, NP] = generateSyntheticData(nmax)
    % Generate synthetic financial data using GARCH models
    % Inputs:
    %   nmax: number of data points to generate
    % Outputs:
    %   X: normalized returns matrix
    %   xav: mean of each asset's returns
    %   NP: number of assets

    % Define GARCH models
    models = {
        garch('Constant',0.01,'GARCH',0.1,'ARCH',0.1)
        garch('Constant',0.01,'GARCH',0.2,'ARCH',0.2)
        garch('Constant',0.01,'GARCH',0.3,'ARCH',0.3)
        garch('Constant',0.01,'GARCH',0.4,'ARCH',0.4)
        garch('Constant',0.01,'GARCH',0.5,'ARCH',0.1)
        garch('Constant',0.01,'GARCH',0.1,'ARCH',0.7)
    };
    
    % Simulate returns
    returns = zeros(nmax, length(models));
    for i = 1:length(models)
        [~, returns(:,i)] = simulate(models{i}, nmax);
    end
    
    % Normalize data
    X = returns;
    NP = size(X,2);
    xav = zeros(1,NP);
    for i = 1:NP
        xav(i) = mean(X(:,i));
        X(:,i) = (X(:,i)-xav(i))./std(X(:,i));
    end
end