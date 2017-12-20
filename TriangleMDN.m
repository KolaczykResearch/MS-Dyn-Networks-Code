%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Name:    TriangleMDN.m        Created: 11/06/16    Revised: 11/20/16
% 
%% Usage:   Performs a recursive, pyramidal series of local (in scale 
%           and position) penalized likelihood optimizations. These are stored in two separate triangular
%           arrays. Each row of the arrays corresponds to substrings (or 'blocks')
%           of data of length/scale 1,2,...,n.  Additionally, a mask is built 
%           along the way to indicate the positions of 'split' decisions at each
%           scale.
%
%% Inputs:  y_id      := one dimensional time series data
%           Data      := length of the time series
%           p         := df of the Chi-square distribution
%           alpha     := controlled type I error rate
%
%% Output:  lambda := empirical penalty
%% Calls:   Only internal Matlab functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t0, t1, mm] = TriangleMDN(y_id, Data, p, alpha) 

[Y, X] = covar_design(Data, p, y_id);
K = size(Data, 2);
[n, ncov] = size(X);% length of the time series

block = repmat(p, K, 1);
block(y_id) = 0;

lam = 1.5 * (K-1) * p * log(n);

% indicates the minimum obs required is 2^4 (This actually should depend on the number of lags used)
Js = 4;            
t0 = [];
t1 = [];
mm = [];

for j = 1:(2^(Js)-1)
    for i = 1: (n-j+1)
        YY = Y(i:(i+j-1));
        pl = 0.5 * sum(YY.^2);
        t0(i, i+j-1) = pl;  
        t1(i, i+j-1) = pl;
        mm(i, i+j-1) = 0;
    end
end

% Iteratively compute penalized likelihoods over intervals of 
% increasing size, storing the various results in the upper triangles
% of the matrices t0, t1, and mm.

for j = (2^Js):n
    for i = 1: (n-j+1)
        XX = X(i:(i+j-1), :);
        YY = Y(i:(i+j-1));
        len = size(YY, 1);
        tseries = XX * spdiags(1./norms(XX)', 0, ncov, ncov);
        
        lambda = compute_lambda(YY, len, p, K-1, alpha); % Empirical lambda

        [~, obj, ~] = g_lasso(tseries, YY, lambda, block, 1.0, 1.0, true);
        
        pl0 = obj;  % None splitted model
        pl1 = t1(i, i+(1:j-1)-1) + t1(i+(1:j-1), i+j-1)' + lam;   % store the local min of the splited model using j data points into pl1

        t0(i, i+j-1) = pl0;
        t1(i, i+j-1) = min(min(pl1), pl0);      % compare the nonsplited model (pl0) versus the smallest among the splitted model (min(pl1))
        [~, mm(i, i+j-1)] = min(pl1);           % mm tracks where the min of pl1 occurs
        mm(i, i+j-1) = mm(i, i+j-1) + i - 1;    % shift the position by i to correct it    
    end
end


    



