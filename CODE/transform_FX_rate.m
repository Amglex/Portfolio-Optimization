function [FX_rate] = transform_FX_rate(FX_rate)

% Transformation of FX-rates
% Function to directly transform data sets that are loaded via the function
% 'hist_stock_data'. This shall provide the most efficient way to proceed 
% with loaded data.
%
% INPUT         FX_rate ... structure
%
% OUTPUT        FX_rate ... matrix with FX rates
%
% MATLAB project, alexander.steeb@student.unisg.ch,
%                 jonas.gartenmeier@student.unisg.ch
% 14.12.2017

 
% Convert structure to table.
data_FX_rate = struct2table(FX_rate);
 
% Vectorize ticker names.
ticker_data_FX_rate = data_FX_rate.Ticker;
 
% Calculate the number of stocks.
n_FX_rate = length(ticker_data_FX_rate);
 
% Retrieve matrix with weekly FX rates depending on number of stocks.
for i = 1:n_FX_rate
     
    FX_rate_weekly(i) = data_FX_rate{i, 'AdjClose'};
     
end
 
% Transform cells to normal matrix.
FX_rate = cell2mat(FX_rate_weekly);
 
end % of function.

