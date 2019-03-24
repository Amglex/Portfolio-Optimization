function [return_ef_classic, std_ef_classic] = calc_ef_classic(mu, sigma)

% Calculation of the efficient frontier [classical way]
% This function calculates the optimal weights for a given mean and sigma 
% for a range of target values. The output consists of two vectors of the
% same length that contain the expected return and standard deviation of
% the different optimal portfolios on the efficient frontier.
%
% INPUT         mu     1xN ... mean-return vector 
%               sigma  NxN ... variance-covariance matrix
%
% OUTPUT        return_ef_classic 1xN ... expected return vector
%               std_ef_classic    1xN ... expected standard deviation
%
% MATLAB project, alexander.steeb@student.unisg.ch,
%                 jonas.gartenmeier@student.unisg.ch
% 14.12.2017

% Number of stocks.
n_stocks = length(mu);

% Vector of suitable target returns based on the average of the expected
% return vector. 
min = mean(mu)* -1; 
max = mean(mu)* 10;

target = linspace(min,max);

% The function calc_opt_weight is used to calculate a matrix consisting of
% the optimal weights for all values of the target return vector.
i = 1;
for t = target
    
    weight_temporary = calc_opt_weight(mu, sigma, t);
    weight_optimal(i,1:n_stocks) = weight_temporary.';
    i = i + 1;
    
end

% For each optimal portfolio the expected return is calculated.
return_ef_classic = weight_optimal * mu.';


% For each optimal portfolio the standard deviation is calcuated by
% multiplying the square root of the optimal portfolio weight with its
% transposed and the respective variance-covariance matrix.
for k = 1:size(weight_optimal,1)
    
    std_ef_classic(k) = sqrt(weight_optimal(k,:) * sigma * weight_optimal(k,:).');
    
end

end % of function

