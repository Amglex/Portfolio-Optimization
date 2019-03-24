function [target_return_tp, tangency_line_tp] = calc_tp(mu,sigma)

% Calculation of tangency portfolio 
% This function calculates the tangency portfolio for a given 'mu'
% and'sigma'.
%
% INPUT         mu     1xN ... mean-return vector 
%               sigma  NxN ... variance-covariance matrix
%
% OUTPUT        target_return_tp 	1xN ... expected return vector
%               tangency_line_tp    	1xN ... expected standard deviation
%
% MATLAB project, alexander.steeb@student.unisg.ch,
%                 jonas.gartenmeier@student.unisg.ch
% 14.12.2017


% A risk free rate of 1% is assumed and calculated on a weekly basis.
risk_free_rate = 0.01/52;

sigma_tp = sigma;
sigma_tp_i = inv(sigma_tp);
mu_tp = mu;
u_tp = ones(1,size(mu_tp,2));

% Vector of suitable target returns based on the average return vector.
min = 0;
max = mean(mu)* 35;

target_return_tp = linspace(min,max);

X = sigma_tp_i * (mu_tp - risk_free_rate .* u_tp ).';

w_tp = X / ( u_tp * X );

return_tp = mu_tp * w_tp;

std_tp = sqrt(w_tp.' * sigma_tp * w_tp);

tangency_line_tp = risk_free_rate + ((return_tp - risk_free_rate ) / std_tp ) * target_return_tp;

end

