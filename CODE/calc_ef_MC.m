function [return_ef_MC, std_ef_MC] = calc_ef_MC(mu, sigma, steps, iterations)

% Calculation of the efficient frontier [Monte Carlo simulation]
% This function calculates the efficient frontier based on the Monte-Carlo 
% method. Unlike the classical efficient frontier this function uses the 
% Monte-Carlo method to calculate the optimal portfolio for each of the 
% target returns. Two additional inputs are therefore needed. The input 
% steps refers to the number of random returns that are generated for each 
% asset. The input iterations refers to how many portfolios are calculated 
% based on those randomly generated returns. 
%
% INPUT         mu              1xN ... mean-return vector 
%               sigma           NxN ... variance-covariance matrix
%               steps           1x1 ... number of steps  
%               iterations      1x1 ... number of iterations
%
% OUTPUT        return_ef_MC    1xN ... expected return vector
%               std_ef_MC       1xN ... expected standard deviation
%
% MATLAB project, alexander.steeb@student.unisg.ch,
%                 jonas.gartenmeier@student.unisg.ch
% 14.12.2017


% Defining the number of paths.
paths = size(mu,2);

% Vector of suitable target returns based on the average of the expected
% return vector. 
min = mean(mu)* -7; 
max = mean(mu)* 15;

target = linspace(min,max);

% In order to receive a matrix of simulated returns (255 * 30), we have to 
% decompose the variance-covariance matrix we obtained previously by using 
% the Cholesky decomposition. 
decomp = chol(sigma,'lower');

% For a pre-defined number of iterations and target, we simulate the
% target returns. 
weight_MC = zeros(size(target,2), paths);

index = 1;

for t = target 

    % Predefine matrix to optimize computation by saving unnecessary steps.
    optimal_weight_MC = zeros(iterations, paths);
    
    % For a pre-defined numbers of iterations and paths we calcuate the 
    % optimal portfolios based on randomly generated returns. 
    for i = 1 : iterations

        random_numbers = randn(steps, paths);

        return_MC = mu + (decomp * random_numbers.').';

        [mu_MC, sigma_MC] = calc_mu_sigma(return_MC);

        weight_temporary_MC = calc_opt_weight(mu_MC, sigma_MC, t);
        
        optimal_weight_MC(i,1:paths) = weight_temporary_MC;

    end

% Save average portfolio weight over all iterations in new matrix.
weight_MC(index,1:paths) = mean(optimal_weight_MC);   

index = index + 1;    

end

% Calculate matrix of returns.
return_ef_MC = weight_MC * mu.';

% Calculate vector of standard deviation.
for k = 1:size(weight_MC,1)
    
    std_ef_MC(k) = sqrt(weight_MC(k,:) * sigma * weight_MC(k,:).');
    
end

end % of function.