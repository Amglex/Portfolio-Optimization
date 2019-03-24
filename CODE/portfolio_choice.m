%% MATLAB Homework - Portfolio Choice, Exercise 8

%  Authors: Alexander Steeb, 17-614-611, alexander.steeb@student.unisg.ch
%           Jonas Gartenmeier, 13-612-700, jonas.gartenmeier@student.unisg.ch
%
%  Date:    14 December 2017
%
%  Outline of different parts
% 
%  1. Loading stock data and FX-rates for calculation of the efficient 
%     frontier (via the classical way and Monte Carlo simulation).
%  2. Set up of US and German portfolio according to pre-specified
%     investment universe.
%  3. Calculation of tangency portfolio according to the last trading days
%     of the years 2013 to 2016.

%% General settings

% Setting random seed to make results reproducable. 
rng(1) 

% Clear all.
clear

%% Declaration of stocks and FX rates ticker symbols

% We included all stocks from the DJIA and the DAX as of 1st of November 
% 2017 via their individual ticker symbols.We excluded, however, one stock 
% from the DAX, namely Vonovia, due to its IPO in July, 2013. We applied
% the same structure to the FX rates. FX rates are EUR to USD, and vice 
% versa, USD to EUR. 
DJIA_stock = {'AAPL';'AXP';'BA';'CAT';'CSCO';'CVX';'KO';'DWDP';'XOM';'GE';'GS';'HD';'IBM';'INTC';'JNJ';'JPM';'MCD';'MMM';'MRK';'MSFT';'PFE';'NKE';'PG';'TRV';'UNH';'UTX';'V';'VZ';'WMT';'DIS'};

DAX_stock = {'CBK.DE';'DBK.DE';'DAI.DE';'ADS.DE';'HEI.DE';'MRK.DE';'SAP.DE';'DPW.DE';'RWE.DE';'TKA.DE';'IFX.DE';'FME.DE';'BAYN.DE';'ALV.DE';'LHA.DE';'DTE.DE';'CON.DE';'FRE.DE';'MUV2.DE';'EOAN.DE';'BAS.DE';'BMW.DE';'SIE.DE';'VOW3.DE';'LIN.DE';'HEN3.DE';'BEI.DE';'DB1.DE';'PSM.DE'};

FX_ticker = {'EURUSD=X';'USDEUR=X'};

%% PART 1
%  ========================================================================

%% Loading the data for the 30 DJIA stocks (daily)

% Retrieve daily data of pre-defined DJIA stocks as structured array from
% Yahoo Finance for a given time period. Code for function
% 'get_yahoo_stockdata3' by Captain Awesome.
DJIA_data_daily = get_yahoo_stockdata3(DJIA_stock,'01-Nov-2012','01-Nov-2017','d');

%% Data transformation

% Transformation of retrieved data into matrix containing only adjusted 
% stock prices.
[stock_price_DJIA] = transform_stock_data(DJIA_data_daily);

%% Calculation of continuous returns

% Calculation of continous returns, therefore using the logarithm.
stock_return_DJIA = log(stock_price_DJIA(2:end,:) ./ stock_price_DJIA(1:end-1,:));

% Plotting continuous returns.
plot(stock_return_DJIA(:,:));
xlabel('Days [01.11.2012 - 01.11.2017]')
ylabel('Daily return')
title('Daily returns for the DJIA-listed stocks');

%% Inspecting the data

% Looking for NaN in the retrieved data set in order to avoid errors in the
% following calculations.
NaN_DJIA_data = find(isnan(stock_return_DJIA));

% If NaN is found, NaN is replaced with 0.
stock_return_DJIA(NaN_DJIA_data) = 0;

%% Calculation of efficient frontier

% Function for obtaining the mean return "mu" and variance-covariance
% matrix "sigma".
[mu, sigma] = calc_mu_sigma(stock_return_DJIA);

% Function for obtaining the efficient frontier the classical way.
[return_ef_classic, std_ef_classic] = calc_ef_classic(mu, sigma);

%% Plotting continuous returns and efficient frontier

% Scatter plot of annualized, continuous returns.
scatter(sqrt(diag(sigma))*sqrt(252), mu*252);
hold on

% Displacement so the text does not overlay the data points.
dx = 0.003;
dy = 0.002;

% Naming data points according to their ticker symbol.
text(sqrt(diag(sigma))*sqrt(252)+dy,mu*252+dx,DJIA_stock);
hold on

% Displaying efficient frontier (annualized, continuous returns) according 
% to industry standard of 252 trading days.
plot(sqrt(252)*(std_ef_classic),return_ef_classic*252,'-')
xlabel('Standard deviation')
ylabel('Annualized return')
title('Scatterplot and efficient frontier of annualized, daily returns of DJIA stocks');

% Limit the range of the x- and y-axis accordingly.
xlim([0.05 0.3])
ylim([-0.1 0.4])

%% Monte Carlo simulation

% Clear certain variables.
clear dx dy NaN_DJIA_data

% Function for obtaining the efficient frontier via Monte Carlo simulation.
[return_ef_MC, std_ef_MC] = calc_ef_MC(mu, sigma, 252, 100);

%% Plotting efficient frontier of Monte Carlo simulation

% Displaying efficient frontier (annualized, continuous returns) according 
% to industry standard of 252 trading days.
plot(sqrt(252)*(std_ef_MC),return_ef_MC*252,'-')
xlabel('Standard deviation')
ylabel('Annualized return')
title('Simulated efficient frontier of annualized, daily returns of DJIA stocks');

%% Random selection of subset of 5 stocks from the DJIA

% Selection of 5 stocks, sampled uniformly at random from the integers 1 to
% 30 (according to the 30 stocks in the DJIA).
r = randsample(29,5);
 
% Ticker symbols of randomly selected stocks.
random_stock = DJIA_stock(r);

% Stock prices of the 5, randomly selected stocks.
stock_random_DJIA = stock_return_DJIA(:,r);

%% Calculation of efficient frontier for randomly selected stocks

% Function for obtaining the mean return "mu" and variance-covariance
% matrix "sigma" for the randomly selected stocks.
[mu_random, sigma_random] = calc_mu_sigma(stock_random_DJIA);

% Function for obtaining the efficient frontier [the classical way] for the
% randomly selected stocks.
[return_random_ef, std_random_ef] = calc_ef_classic(mu_random, sigma_random);

%% Plotting efficient frontier [the classical way] for randomly selected stocks

% Display efficient frontier (annualized, continuous returns) according to 
% industry standard of 252 trading days.
plot(sqrt(252)*(std_random_ef),return_random_ef*252,'-')
xlabel('Standard deviation')
ylabel('Annualized return')
title('Efficient frontier of annualized, daily returns for a subset of 5 DJIA stocks');

%% Monte Carlo simulation for subset of 5 randomly selected stocks 

% Function for obtaining the efficient frontier of the 5 randomly selected 
% stock via Monte Carlo simulation.
[return_random_MC, std_random_MC] = calc_ef_MC(mu_random, sigma_random, 252, 100);

%% Plotting efficient frontier of Monte Carlo simulation for subset of 5 randomly selected stocks

% Display efficient frontier (annualized, continuous returns) according to 
% industry standard of 252 trading days.
plot(sqrt(252)*(std_random_MC),return_random_MC*252,'-')
xlabel('Standard deviation')
ylabel('Annualized return')
title('Simulated efficient frontier of annualized, daily returns for a subset of 5 DJIA stocks');

%% PART 2
%  ========================================================================

%% General setting

% Clear all variables except given tickers.
clearvars -except DAX_stock DJIA_stock FX_ticker

%% Loading and transforming the data of the 30 DJIA stocks (weekly)

% Retrieve weekly data of pre-defined DJIA stocks as structured array from
% Yahoo Finance for a given time period.
DJIA_data_weekly = get_yahoo_stockdata3(DJIA_stock,'01-Nov-2012','01-Nov-2017','w');

% Transformation of retrieved data into matrix containing only adjusted 
% stock prices.
[stock_price_DJIA_weekly] = transform_stock_data(DJIA_data_weekly);

%% Loading and transforming the data of the 29 DAX stocks (weekly)

% Retrieve weekly data of pre-defined DAX stocks as structured array from
% Yahoo Finance for a given time period.
DAX_data_weekly = get_yahoo_stockdata3(DAX_stock,'01-Nov-2012','01-Nov-2017','w');

% Transformation of retrieved data into matrix containing only adjusted 
% stock prices.
[stock_price_DAX_weekly] = transform_stock_data(DAX_data_weekly);

%% Loading and transforming the data of the FX rates (weekly)

% Retrieve weekly data of exchange rates as structured array from
% Yahoo Finance for a given time period. Code for function
% 'hist_stock_data' by Josiah Renfree.
FX_data_weekly = hist_stock_data('01112012','01112017', FX_ticker,'frequency','wk');

% Transformation of retrieved data into matrix containing only adjusted 
% stock prices.
[FX_rate_weekly] = transform_FX_rate(FX_data_weekly);

%% Calculation of continuous returns for US portfolio

% Conversion of weekly stock prices of each DAX stock from EUR to USD.
stock_price_DAX_USD = FX_rate_weekly(:,1) .* stock_price_DAX_weekly;

% Set up of US investment portfolio [30 DJIA stocks, 29 DAX stocks, denoted
% in USD].
US_port = [stock_price_DJIA_weekly, stock_price_DAX_USD];
    
% Calculating continuous returns of portfolio held by US investor.
return_US_port = log(US_port(2:end,:) ./ US_port(1:end-1,:));

% Plotting continuous returns.
plot(return_US_port(:,:))
xlabel('Week')
ylabel('Weekly return')
title('Weekly returns for US portfolio [262 weeks]');

%% Inspecting the data of US portfolio

% Looking for NaN in the retrieved data set in order to avoid errors in 
% the following calculations.
NaN_US_port = find(isnan(return_US_port));

% If NaN is found, NaN is replaced with 0.
return_US_port(NaN_US_port) = 0;

%% Calculation of continuous returns for GER portfolio

% Conversion of weekly stock prices of each DJIA stock from USD to EUR.
stock_price_DJIA_EUR = FX_rate_weekly(:,2) .* stock_price_DJIA_weekly;

% Set up of GER investment portfolio [30 DJIA stocks, 29 DAX stocks, denoted
% in EUR].
GER_port = [stock_price_DAX_weekly, stock_price_DJIA_EUR];

% Calculating continuous returns of portfolio held by GER investor.
return_GER_port = log(GER_port(2:end,:) ./ GER_port(1:end-1,:));

% Plotting continuous returns.
plot(return_GER_port(:,:))
xlabel('Week')
ylabel('Weekly return')
title('Weekly returns for GER portfolio [262 weeks]');

%% Inspecting the data of GER portfolio

% Looking for NaN in the retrieved data set in order to avoid errors in 
% the following calculations.
NaN_GER_port = find(isnan(return_GER_port));

% If NaN is found, NaN is replaced with 0.
return_GER_port(NaN_GER_port) = 0;

%% Calculation of efficient frontier for US portfolio and GER portfolio

% Function for obtaining the mean return "mu" and variance-covariance
% matrix "sigma" for the US investor.
[mu_US_port, sigma_US_port] = calc_mu_sigma(return_US_port);

% Function for obtaining the Efficient Frontier [the classical way] for the
% US investor.
[return_US_port_ef, std_US_port_ef] = calc_ef_classic(mu_US_port, sigma_US_port);

%% German Investor

% Function for obtaining the mean return "mu" and variance-covariance
% matrix "sigma" for German investor.
[mu_GER_port, sigma_GER_port] = calc_mu_sigma(return_GER_port);

% Function for obtaining the efficient frontier [the classical way] for 
% German investor.
[return_GER_port_ef, std_GER_port_ef] = calc_ef_classic(mu_GER_port, sigma_GER_port);

% Displaying efficient frontier (annualized, continuous returns) for German
% and US portfolio.
plot(sqrt(12)*(std_GER_port_ef),return_GER_port_ef*12,'-b', sqrt(12)*(std_US_port_ef),return_US_port_ef*12,'-r')
xlabel('Standard deviation')
ylabel('Return')
legend('German Investor','US Investor')
title('Efficient frontier of annualized, weekly returns for GER portfolio');

%% %% PART 3
%  ========================================================================

%% General setting

% Clear all variables except given tickers.
clearvars -except DAX_stock DJIA_stock FX_ticker

%% Loading the data for the tangency portfolios

% In order to compute the tangency portfolio for the last trading day of 
% the years 2016, 2015 and 2014 (each with a time period of 3 years) weekly
% rates of the DJIA for the individual time period have to be loaded. 
DJIA_data_tp_1 = get_yahoo_stockdata3(DJIA_stock,'01-Jan-2013','31-Dez-2016','w');

DJIA_data_tp_2 = get_yahoo_stockdata3(DJIA_stock,'01-Jan-2012','31-Dez-2015','w');

DJIA_data_tp_3 = get_yahoo_stockdata3(DJIA_stock,'01-Jan-2011','31-Dez-2014','w');

%% Transforming the data sets and computing the continuous returns for each tangency portfolio.

% Transform each data set in order to retrieve adjusted stock prices.
[stock_price_DJIA_tp_1] = transform_stock_data(DJIA_data_tp_1);

[stock_price_DJIA_tp_2] = transform_stock_data(DJIA_data_tp_2);

[stock_price_DJIA_tp_3] = transform_stock_data(DJIA_data_tp_3);

% Compute the stock returns for each data set.
stock_return_DJIA_tp_1 = log(stock_price_DJIA_tp_1(2:end,:) ./ stock_price_DJIA_tp_1(1:end-1,:));

stock_return_DJIA_tp_2 = log(stock_price_DJIA_tp_2(2:end,:) ./ stock_price_DJIA_tp_2(1:end-1,:));

stock_return_DJIA_tp_3 = log(stock_price_DJIA_tp_3(2:end,:) ./ stock_price_DJIA_tp_3(1:end-1,:));

%% Inspecting the data of tangency portfolios

% Looking for NaN in the retrieved data set in order to avoid errors in 
% the following calculations.
NaN_DJIA_data_tp_1 = find(isnan(stock_return_DJIA_tp_1));

NaN_DJIA_data_tp_2 = find(isnan(stock_return_DJIA_tp_2));

NaN_DJIA_data_tp_3 = find(isnan(stock_return_DJIA_tp_3));

% If NaN is found, NaN is replaced with 0.
stock_return_DJIA_tp_1(NaN_DJIA_data_tp_1) = 0;

stock_return_DJIA_tp_2(NaN_DJIA_data_tp_2) = 0;

stock_return_DJIA_tp_3(NaN_DJIA_data_tp_3) = 0;

%% Calculation of efficient frontier for each tangency portfolio

% Function for obtaining the mean return "mu" and variance-covariance
% matrix "sigma" for each tangency portfolio.
[mu_tp_1, sigma_tp_1] = calc_mu_sigma(stock_return_DJIA_tp_1);

[mu_tp_2, sigma_tp_2] = calc_mu_sigma(stock_return_DJIA_tp_2);

[mu_tp_3, sigma_tp_3] = calc_mu_sigma(stock_return_DJIA_tp_3);

% Function for obtaining the efficient frontier [the classical] way for
% each tangency portfolio.
[return_tp_1, std_tp_1] = calc_ef_classic(mu_tp_1, sigma_tp_1);

[return_tp_2, std_tp_2] = calc_ef_classic(mu_tp_2, sigma_tp_2);

[return_tp_3, std_tp_3] = calc_ef_classic(mu_tp_3, sigma_tp_3);

% Calculation and plotting of tangency line for each tangency
% portfolio with its respective efficient frontier.
[target_return_tp_1, tangency_line_tp_1] = calc_tp(mu_tp_1,sigma_tp_1);

[target_return_tp_2, tangency_line_tp_2] = calc_tp(mu_tp_2,sigma_tp_2);

[target_return_tp_3, tangency_line_tp_3] = calc_tp(mu_tp_3,sigma_tp_3);

plot(sqrt(52)*target_return_tp_1, tangency_line_tp_1*52,'-r',sqrt(52)*(std_tp_1),return_tp_1*52,'-r')
hold on
plot(sqrt(52)*target_return_tp_2, tangency_line_tp_2*52,'-g',sqrt(52)*(std_tp_2),return_tp_2*52,'-g')
hold on
plot(sqrt(52)*target_return_tp_3, tangency_line_tp_3*52,'-b',sqrt(52)*(std_tp_3),return_tp_3*52,'-b')
xlabel('Standard Deviation')
ylabel('Return')
legend('Tangency Portfolio 01.01.13 - 31.12.16','Tangency Portfolio 01.01.12 - 31.12.15','Tangency Portfolio 01.01.11 - 31.12.14')
title('Tangency portfolios over time')

%% Declaration

% We hereby certify that
% We have written the program ourselves except for clearly marked pieces of
% code. We have tested the program and it ran without crashing (if applicable)
% (Alexander Steeb, Jonas Gartenmeier)
