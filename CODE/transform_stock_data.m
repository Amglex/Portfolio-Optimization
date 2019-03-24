function [stock_prices] = transform_stock_data(stock_data)

% Transformation of stock data
% Function to directly transform data sets that are loaded via the function
% 'get_yahoo_stockdata3'. This shall provide the most efficient way to proceed 
% with loaded data.
%
% INPUT         stock_data     	... structured array 
%
% OUTPUT        stock_prices    NxN ... matrix with stock prices
%
% MATLAB project, alexander.steeb@student.unisg.ch,
%                 jonas.gartenmeier@student.unisg.ch
% 14.12.2017



% Convert structure to matrix according to the length of loaded data.
for i = 1:length(stock_data)
    
    % Retrieve the adjusted closing prices instead of closing prices in
    % order to avoid unexpected price fluctuations due to stock splits,
    % etc.
    stock_prices(:,i) = stock_data{i,1}.adjClosePrice;
    
end % of function

