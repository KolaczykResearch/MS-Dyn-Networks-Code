%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Name:    compute_lambda.m        Created: 11/06/16    Revised: 11/20/16
% 
%% Usage:   Compute the empirical penalty for the group lasso. 
%           'alpha' controls the Type I error of recovering the connected components
%
%% Inputs:  Y      := one dimensional time series data
%           T      := length of the time series
%           p      := df of the Chi-square distribution
%           N      := number of nodes in the neighborhood
%           alpha  := controlled type I error rate
%
%% Output:  lambda := empirical penalty
%% Calls:   Only internal Matlab functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lambda] = compute_lambda(Y, T, p, N, alpha)

    lambda =  sqrt(sum(Y.^2)/T) * sqrt(p * chi2inv(1 - alpha/(2 * N*(N-1)), p));
    
end
