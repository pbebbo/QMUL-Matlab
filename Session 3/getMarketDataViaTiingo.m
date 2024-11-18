function [data, price_series, dates] = getMarketDataViaTiingo(symbols, startdate, enddate, interval)
    % Downloads market data from Tiingo for multiple ticker symbols
    % 
    % INPUT:
    % symbols   - cell array of ticker symbols i.e. {'AAPL', 'MSFT'}
    % startdate - start date for data request
    % enddate   - end date for data request
    % interval  - data frequency (daily, weekly, monthly)
    %
    % OUTPUT:
    % data         - struct array with data for each symbol
    % price_series - matrix of adjusted close prices (columns are symbols)
    % dates        - matrix of dates in datenum format
    
    % Handle input parameters
    if ~iscell(symbols)
        symbols = {symbols};  % Convert single string to cell array
    end
    
    if(nargin() == 1)
        startdate = datetime('1-Jan-2018');
        enddate = datetime('today');
        interval = 'daily';
    elseif (nargin() == 2)
        enddate = datetime('today');
        interval = 'daily';
    elseif (nargin() == 3)
        interval = 'daily';
    end
    
    % Initialize outputs
    data(length(symbols)) = struct('Ticker', [], 'Date', [], 'Open', [], ...
                                  'High', [], 'Low', [], 'Close', [], ...
                                  'Volume', [], 'AdjClose', []);
    price_series = [];
    dates = [];
    
    % Get API key from environment variable
    api_key = getenv('TIINGO_API_KEY');
    if isempty(api_key)
        error('TIINGO_API_KEY environment variable not set');
    end
    
    % Format dates for API
    startStr = datestr(startdate, 'yyyy-mm-dd');
    endStr = datestr(enddate, 'yyyy-mm-dd');
    
    % Process each symbol
    for i = 1:length(symbols)
        symbol = symbols{i};
        
        % Construct URL
        baseUrl = 'https://api.tiingo.com/tiingo/daily/%s/prices';
        url = sprintf([baseUrl '?startDate=%s&endDate=%s&resampleFreq=%s'], ...
                     upper(symbol), startStr, endStr, interval);
        
        % Set up options
        options = weboptions('ContentType', 'json', ...
                            'HeaderFields', {...
                                'Authorization', ['Token ' api_key]});
        
        try
            % Download data
            response = webread(url, options);
            
            if ~isempty(response)
                % Extract dates
                dates_dt = datetime({response.date}, ...
                           'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS''Z');
                
                % Fill the struct fields
                data(i).Ticker = symbol;
                data(i).Date = datenum(dates_dt);
                data(i).Open = [response.open]';
                data(i).High = [response.high]';
                data(i).Low = [response.low]';
                data(i).Close = [response.close]';
                data(i).Volume = [response.volume]';
                data(i).AdjClose = [response.adjClose]';
                
                % Add to price series and dates matrices
                price_series(:,i) = data(i).AdjClose;
                dates(:,i) = data(i).Date;
            else
                warning('No data available for %s', symbol);
            end
            
        catch ME
            warning('Failed to download data for %s: %s', symbol, ME.message);
        end
    end
end