%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Name:    RecoverCoef.m        Created: 11/06/16    Revised: 11/20/16
% 
%% Usage:   recover the piecewise stationary AR coefficients 
%
%% Inputs:  Data   := multi-dimensional time series data
%           bp     := detected change points
%           alpha  := controlled type I error rate
%           p      := number of lags
%
%% Output:  coef := coefficient
%% Calls:   Only internal Matlab functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [coef] = RecoverCoef(Data, bp, alpha, p)

coef = [];
[n, K] = size(Data); 

bp = [1 bp n-p];
num_intvl = length(bp)-1;

for k = 1:num_intvl,            % loop throught each stationary blocks
    nrep = bp(k+1)-bp(k);
    xind = bp(k):bp(k+1);
    oneslide = [];
    for y_id = 1:K,
        tmpcoef = [];
        block = repmat(p, K, 1);
        block(y_id) = 0;
        [Y, X] = covar_design(Data, p, y_id);
        [~, ncov] = size(X);  

        tseries = X(xind, :);
        Ylocal = Y(xind, :);
        n = size(Ylocal, 1);
        tseries = tseries * spdiags(1./norms(tseries)', 0, ncov, ncov); % normalize covariates

        lambda = compute_lambda(Ylocal, n, p, K-1, alpha)*1/k^1.5;      % compute empirical lambda 

        [x, ~, ~] = g_lasso(tseries, Ylocal, lambda, block, 1.0, 1.0, true);
        x = x(1:p:(K-1)*p);
        
        tmpcoef = [tmpcoef x];
        tmpcoef = insertrows(tmpcoef, zeros(1, 1), y_id-1);
        oneslide = [oneslide tmpcoef];
    end
    coef = cat(3, coef, reshape(repmat(oneslide, 1, nrep), K, K, nrep)); % reshape the coefficient matrix
end