function stock_data = getMarketDataViaTiingo(symbols, startDate, endDate, freq)
% getMarketDataViaTiingo  Fetch end-of-day prices from Tiingo
%
% stock_data = getMarketDataViaTiingo(symbols, startDate, endDate, freq)
% Returns struct array with fields:
%   .ticker (string)
%   .data   (timetable with variables open, high, low, close, volume,
%           adjOpen, adjHigh, adjLow, adjClose, adjVolume, divCash, splitFactor)
arguments
    symbols {mustBeNonempty}
    startDate
    endDate
    freq (1,:) char {mustBeMember(freq,{'daily'})} = 'daily'
end

% Normalize inputs
if ischar(symbols) || isstring(symbols)
    symbols = cellstr(symbols);
end

startDT = datetime(startDate);
endDT   = datetime(endDate);

apiKey = getenv('TIINGO_API_KEY');
if isempty(apiKey)
    error('TIINGO_API_KEY environment variable not set. Use setenv(''TIINGO_API_KEY'',''<key>'').');
end

baseUrl = 'https://api.tiingo.com/tiingo/daily/';
opts = weboptions('Timeout', 30, ...
                  'HeaderFields', {'Authorization', ['Token ' apiKey]}, ...
                  'ContentType', 'json');

n = numel(symbols);
stock_data = repmat(struct('ticker', "", 'data', timetable.empty), n, 1);

for i = 1:n
    t = strtrim(symbols{i});
    stock_data(i).ticker = string(t);
    
    try
        url = sprintf('%s%s/prices?startDate=%s&endDate=%s', ...
            baseUrl, upper(t), ...
            datestr(startDT, 'yyyy-mm-dd'), ...
            datestr(endDT, 'yyyy-mm-dd'));
        
        fprintf('Fetching %s... ', t);
        resp = webread(url, opts);
        
        if isempty(resp)
            fprintf('No data returned\n');
            continue
        end
        
        % Convert struct array to timetable
        nRows = numel(resp);
        fprintf('%d rows retrieved... ', nRows);
        
        % Pre-allocate arrays
        dates = NaT(nRows, 1);
        open = nan(nRows, 1);
        high = nan(nRows, 1);
        low = nan(nRows, 1);
        closep = nan(nRows, 1);
        volume = nan(nRows, 1);
        adjOpen = nan(nRows, 1);
        adjHigh = nan(nRows, 1);
        adjLow = nan(nRows, 1);
        adjClose = nan(nRows, 1);
        adjVolume = nan(nRows, 1);
        divCash = nan(nRows, 1);
        splitFactor = nan(nRows, 1);
        
        % Extract data row by row
        for j = 1:nRows
            % Parse date and remove timezone to avoid conflicts
            dt = datetime(resp(j).date, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS''Z''');
            dates(j) = datetime(dt, 'TimeZone', '');  % Remove timezone
            
            open(j) = safeExtract(resp(j), 'open');
            high(j) = safeExtract(resp(j), 'high');
            low(j) = safeExtract(resp(j), 'low');
            closep(j) = safeExtract(resp(j), 'close');
            volume(j) = safeExtract(resp(j), 'volume');
            adjOpen(j) = safeExtract(resp(j), 'adjOpen');
            adjHigh(j) = safeExtract(resp(j), 'adjHigh');
            adjLow(j) = safeExtract(resp(j), 'adjLow');
            adjClose(j) = safeExtract(resp(j), 'adjClose');
            adjVolume(j) = safeExtract(resp(j), 'adjVolume');
            divCash(j) = safeExtract(resp(j), 'divCash');
            splitFactor(j) = safeExtract(resp(j), 'splitFactor');
        end
        
        % Build timetable
        T = table(open, high, low, closep, volume, ...
                  adjOpen, adjHigh, adjLow, adjClose, adjVolume, ...
                  divCash, splitFactor, ...
                  'VariableNames', {'open', 'high', 'low', 'close', 'volume', ...
                                   'adjOpen', 'adjHigh', 'adjLow', 'adjClose', 'adjVolume', ...
                                   'divCash', 'splitFactor'});
        
        tt = table2timetable(T, 'RowTimes', dates);
        stock_data(i).data = sortrows(tt);
        
        fprintf('SUCCESS\n');
        
    catch ME
        fprintf('FAILED: %s\n', ME.message);
        warning('Failed to fetch %s: %s', t, ME.message);
    end
end
end

function val = safeExtract(s, fieldName)
% Safely extract numeric value from struct field
    if isfield(s, fieldName) && ~isempty(s.(fieldName))
        v = s.(fieldName);
        if isnumeric(v)
            val = double(v);
        else
            val = NaN;
        end
    else
        val = NaN;
    end
end