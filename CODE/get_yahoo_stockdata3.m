function stock = get_yahoo_stockdata3(ticker,d1,d2,freq)
% Code by Captain Awesome
%
% Updated from v3 when in May 2017, yahoo went and changed how stock data
% was shown on web pages. 
%
% INPUTS:
%
% ticker  <-- Yahoo ticker symbol for desired security. This can be a char
%             string for a single stock, or data can be retrieved for
%             multiple stocks by using a cellstr array.
%
% d1      <-- start date for data. Can be a matlab datenumber or a date string.
%             Default = 100 days ago
%
% d2      <-- end date for data. Can be a matlab datenumber or a date string.
%             Default = today
%
% freq    <-- data frequency 'd' (daily), 'w' (weekly), or 'm' (monthly).
%             Default = 'd'
%
% OUTPUT:
%
% stock  <-- matlab data structure with output data.
%
% Examples:
%
%  Get data for the past 100 days.
%  stock = get_yahoo_stockdata3('goog');
%  stock = get_yahoo_stockdata3({'goog', 'aapl', 'fb'});
%
%  Get data from 01-Mar-2008 to now.
%  stock = get_yahoo_stockdata3('goog','01-Mar-2008');
%
%  Get data for the past 5 years.
%  stock = get_yahoo_stockdata3('goog', now-5*365, now);
%
%  Get data for specific date range, but weekly instead of daily
%  stock = get_yahoo_stockdata3({'goog', 'aapl', 'fb'},'01-Jan-2009','01-Apr-2010','w');
%
% Captain Awesome, November 2017

if nargin<4
  freq = 'd';
end
if nargin<3
  d2 = now;
end
if nargin<1
  d1 = d2-100;
end
  
d1 = floor(datenum(d1));
d2 = floor(datenum(d2));
ticker = upper(ticker);

if d1>d2
  error('bad date order');
end

if isempty(ticker)
  error('No ticker given.');
end

if sum(strcmpi(freq,{'daily','day','d'}))
  freq = 'd';
elseif sum(strcmpi(freq,{'weekly','week','w','wk'}))
  freq = 'wk';
elseif sum(strcmpi(freq,{'monthly','month','mmowk'}))
  freq = 'mo';  
else
  error('data frequency not recognized');
end

% If given a cellstr array of tickers, then this will recursively call this
% function for each ticker, output will be a cell array of stock data
% structures.
if iscell(ticker)
  stock = cellfun(@(x) get_yahoo_stockdata3(x,d1,d2,freq),...
    ticker, 'uniformoutput', false);
  return
end

clear stockData
stock.ticker     = ticker;
stock.dataSource = 'Yahoo Finance';
stock.dataUpdate = datestr(now,0);
stock.errorMsg   = '';
stock.dataFreq   = freq;

% Yahoo finance uses a unix serial date number, so will have to convert to
% that.  That's a UNIX timestamp -- the number of seconds since January 1, 1970.
unix_epoch = datenum(1970,1,1,0,0,0);
d1u = floor(num2str((datenum(d1) - unix_epoch)*86400));
d2u = floor(num2str((datenum(d2) - unix_epoch)*86400));

site=strcat('https://finance.yahoo.com/quote/',ticker,'/history?',...
  'period1=',d1u,'&period2=',d2u,'&interval=1',freq,'&filter=history&',...
  'frequency=1',freq);

[temp, status] = urlread(site);

if ~status
  warning(['stock data download failed: ',ticker]);
  stock.errorMsg=['stock data download failed (',...
    datestr(now,0),'): ',ticker];
  return
end

%% Check that this is the right ticker data
C = strsplit(temp,'Ta(start) ')';

for k = 1:length(C)
  
  s = C{k};
  
  if sum(strfind(lower(s),lower('ticker=')))
    t = strsplit(s,'ticker=');
    t = t{2};
    t = textscan(t,'%s');
    t = t{1}{1};
    if ~strcmp(t, ticker)
      error('ticker mismatch');
    end
    break
  end
  clear s
  
end

clear k C

%% Get data from 'Historical Prices' section

C = strsplit(temp,'"HistoricalPriceStore":{"prices":[');

if length(C)==1
  % In this case all the data was displayed to the screen and there was no
  % extra data in "HistoricalPriceStore"
  
  data = [];
  ddata = [];
  sdata = [];
  
else
  % In this case only some data was initially displayed and the rest was
  % put in this "HistoricalPriceStore" section and a script would display
  % it as the user scrolled down.
  
  C = strsplit(C{2},'},');
  
  n = length(C);
  data = NaN(n,7);  % stock data
  ddata = NaN(n,3); % dividend data
  sdata = NaN(n,3); % split events
  
  for k=1:n
    
    s = C{k};
    
    if length(s)<13
      continue
    end  
    
    if strcmp(s(1:13),'"eventsData":')
      break    
    end
    
    if sum(strfind(lower(s),lower('"splitRatio"')))

      x = textscan(s,['{"date":%f "numerator":%f "denominator":%f%*s'],...
        'delimiter',',','TreatAsEmpty','null');
    
      if sum(cellfun('isempty',x))
        error('badness');
      end
      
      sdata(k,:) = cell2mat(x);
      clear x   
    
    
    elseif strcmp(s(1:8),'{"date":')
      x = textscan(s,['{"date":%f "open":%f "high":%f "low":%f "close":%f "volume":%f "adjclose":%f}'],...
        'delimiter',',','TreatAsEmpty','null');
      
      if sum(cellfun('isempty',x))
        error('badness');
      end
      
      data(k,:) = cell2mat(x);
      clear x
    
    elseif strcmp(s(1:10),'{"amount":')
      x = textscan(s,['{"amount":%f,"date":%f,"type":"DIVIDEND","data":%f'],...
        'delimiter','','TreatAsEmpty','null');
      
      if sum(cellfun('isempty',x))
        error('badness');
      end
      
      ddata(k,:) = cell2mat(x);
      clear x
      
    end
    
    clear s
    
  end
  clear n k
  
  % data columns: date, open, high, low, close, volume, adjclose
  data(isnan(data(:,1)),:) = [];
  data(:,1) = datenum(datetime(data(:,1),'ConvertFrom','posixtime'));
  data = sortrows(data,1);
  
  % ddate columns: Amount, date, data
  ddata(isnan(ddata(:,1)),:) = [];
  ddata(:,2) = datenum(datetime(ddata(:,2),'ConvertFrom','posixtime'));
  ddata = sortrows(ddata,2);
  
  %sdata columns: Date numerator demonitor
  sdata(isnan(sdata(:,1)),:) = [];
  sdata(:,1) = datenum(datetime(sdata(:,1),'ConvertFrom','posixtime'));
  sdata = sortrows(sdata,2);
  
end

clear C

%% Assign data to output structure

if isempty(ddata)
  stock.dividends = [];
else
  stock.dividends.DateTime = ddata(:,2);
  stock.dividends.Amount = ddata(:,1);
end

if isempty(sdata)
  stock.splits = [];
else
  stock.splits.DateTime = sdata(:,1);
  stock.splits.numerator = sdata(:,2);
  stock.splits.denominator = sdata(:,3);
end


if isempty(data)
  stock.errorMsg=['No data found in stock data download (',datestr(now,0),'): ',ticker];
  warning(['No data found in stock data download: ',ticker]);
  return
end

stock.range = [datestr(data(1,1),1),...
    ' to ',datestr(data(end,1),1)];

stock.varnotes={...
% Variable          Units   Description                     Format 
  'DateTime',      '[EST]', 'Date of stock quote',          'yyyy-mm-dd';...
  'openPrice',     '[$]',   'Opening price of stock',          '%.2f';...
  'highPrice',     '[$]',   'High price of stock',             '%.2f';...
  'lowPrice',      '[$]',   'Low price of stock',              '%.2f';...
  'closePrice',    '[$]',   'Closing price of stock',          '%.2f';...
  'adjClosePrice', '[$]',   'Adjusted close price of stock',   '%.2f';...
  'volume',        '[-]',   'Trading volume',                  '%.0f'};

stock.DateTime      = data(:,1);
stock.openPrice     = data(:,2);
stock.highPrice     = data(:,3);
stock.lowPrice      = data(:,4);
stock.closePrice    = data(:,5);
stock.volume        = data(:,6);
stock.adjClosePrice = data(:,7);


end % function get_yahoo_stockdata3