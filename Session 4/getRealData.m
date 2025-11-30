function [X, xav, Nassets] = getRealData(symbols, startDate, endDate)
    % Fetch and process real market data
    % Inputs:
    %   symbols: cell array of stock symbols
    %   startDate: start date for data fetch
    %   endDate: end date for data fetch
    % Outputs:
    %   X: normalized returns matrix
    %   xav: mean of each asset's returns
    %   Nassets: number of assets

    stock_data = getMarketDataViaTiingo(symbols, startDate, endDate, 'daily');
    
    Ndays = length(stock_data(1).data.adjClose);
    Nassets = length(stock_data);
    X = zeros(Ndays-1, Nassets);
    
    % Calculate returns
    for i = 1:Nassets
        X(:,i) = price2ret(stock_data(i).data.adjClose);
    end
    
    % Normalize data
    xav = zeros(1,Nassets);
    for i = 1:Nassets
        xav(i) = mean(X(:,i));
        X(:,i) = (X(:,i)-xav(i))./std(X(:,i));
    end
end 