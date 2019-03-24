function [mu, sigma] = calc_mu_sigma(stock_return)

% Calculation of mean and sigma.
% This function is used to calculate the mean return vector and the
% variance-covariance matrix for a given matrix of stock returns. 
%
% INPUT         stock_return     ... structured array
%
% OUTPUT        mu           1xN ... mean-return vector
%               sigma        NxN ... variance-covariance matrix
%
% MATLAB project, alexander.steeb@student.unisg.ch,
%                 jonas.gartenmeier@student.unisg.ch
% 14.12.2017

% The two generated outputs are the mean return vector 'mu' and the
% variance-covariance matrix 'sigma'.
mu = mean(stock_return);

sigma = cov(stock_return);

end % of function.

