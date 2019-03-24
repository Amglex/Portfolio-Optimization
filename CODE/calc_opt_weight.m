function [weight_optimal] = calc_opt_weight(mu, sigma, target)

% Calculation of optimal portfolio weights
% This function calculates the optimal portfolio weights for a given mean
% return vector 'mu', the variance-covariance matrix 'sigma' and a given
% target return vector. 
%
% INPUT         mu                1xN ... mean-return vector 
%               sigma             NxN ... variance-covariance matrix
%               t                 1xN ... target return vector
%
% OUTPUT        weight_optimal    1xN ... vector of optimal weights
%
% MATLAB project, alexander.steeb@student.unisg.ch,
%                 jonas.gartenmeier@student.unisg.ch
% 14.12.2017

% The according weight constraint.
u = ones(size(mu));

% Obtaining the inverse of the variance-covariance matrix.
sigma_i = inv(sigma);

% Through constrained optimization with Lagrange multipliers we can
% transofrm the given minimization problem.
A = u * sigma_i * mu.';
B = mu * sigma_i * mu.';
C = u * sigma_i * u.';
D = B * C - A^2;

% Output are the optimal weights.
weight_optimal = 1 / D * (B * sigma_i * u.' - A * sigma_i * mu.' + target .* (C * sigma_i * mu.' - A * sigma_i * u.'));

end % of function.

