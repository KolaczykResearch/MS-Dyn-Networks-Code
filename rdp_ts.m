    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Name:    rdp_ts.m        Created: 11/06/16    Revised: 11/22/16
%
%% Usage:   Fitting the model using recursive dyadic partition.
%
%% Inputs:  Data     := input multivariate time series data
%           p        := number of lags assumed
%           alpha    := type I error rate, used to compute lambda
%
%% Output:  bestFit  := Fitted value
%           decor    := a 0/1 valued sequence indicates the structure of the
%                       model
%           lambdas  := lambdas used in RDP
%           coef     := coefficient used in RDP
%           objr     := value of objection functions in RDP
%% Calls:   Only internal Matlab functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [bestFit, decor, lambdas_whole, coef_whole, objr_whole] = rdp_ts(y_id, Data, p, alpha) 

    [~, K] = size(Data);  % K-1: size of the neighborhood, 
    Js = 4;   
    lambdas_whole = [];
    objr_whole = [];
    coef_whole = [];
    decor = [];
    
    objr = [];
    coef = [];
    lambdas = [];

    [Y, X] = covar_design(Data, p, y_id);
    [n, ncov] = size(X);          % n: length of the time series; ncov: dimension of X
    J = floor(log2(n));

    lam = 0.5*(K-1) * p * log(n);
    bestPL = zeros([2^(J-Js) 1]);
    bestFit = zeros(size(Y));

    block = repmat(p, K, 1);
    block(y_id) = 0;
    %  Starting from the finest interval (with 8 observations)
    for k = 0:1:2^(J-Js) - 1,

        xind = 2^Js * k + 1:2^Js * (k+1);
        tseries = X(xind, :);
        Ylocal = Y(xind, :);

        tseries = tseries * spdiags(1./norms(tseries)', 0, ncov, ncov); % normalize covariates

        lambda = compute_lambda(Ylocal, n, p, K-1, alpha);    % compute empirical lambda 

        [x, obj, ~] = g_lasso(tseries, Ylocal, lambda, block, 1.0, 1.0, true);
        fit = tseries * x;  % compute fitted value

        lambdas = [lambdas lambda];
        bestFit(xind, :) = fit;
        bestPL(k+1) = obj; 
        coef = [coef x];
        objr = [objr obj]; 
    end
    
    objr_whole = [99 objr];
    coef_whole = [99*ones([ncov 1]) coef];

    % Now proceed down through coarser scales.
    for i = J - Js - 1:-1:0, % J-i becomes Js
        decorj = [];
        objr = [];
        coef = [];
        lambdas = [];
        
        oldbestPL = bestPL;
        bestPL = zeros([2^(i) 1]);

        for k = 0:1:2^i-1,
            xind = 2^(J-i)*k+1:2^(J-i)*(k+1);
            plind = [2 * k + 1 2 * k + 2];
            tseries = X(xind,:);
            Ylocal = Y(xind, :);

            len = size(Ylocal, 1);

            tseries = tseries * spdiags(1./norms(tseries)', 0, ncov, ncov); % normalize covariates

            lambda = compute_lambda(Ylocal, len, p, K-1, alpha);    % compute empirical lambda 

            [x, obj, ~] = g_lasso(tseries, Ylocal, lambda, block, 1.0, 1.0, true);
            fit = tseries * x;

            pl0 = obj; % None splitted model
            plsplit = sum(oldbestPL(plind, :)) + lam; % Splitted model

            tmpdecor = (plsplit < pl0); % tmpdecor = 1 means prefer the splitted model, thus split = 1
    %        bestFit(xind, :) = bsxfun(@times, bestFit(xind, :), tmpdecor) + bsxfun(@times, fit, 1-tmpdecor);       
    %        bestPL(k+1,:) = plsplit .* tmpdecor + pl0 .* (1-tmpdecor);
            bestFit(xind) = bestFit(xind) .* tmpdecor + (fit) .* (1 - tmpdecor);

            lambdas = [lambdas lambda];
            bestPL(k+1,:) = plsplit .* tmpdecor + pl0 .* (1-tmpdecor);
            coef = [coef x];
            objr = [objr obj];
            decorj = [decorj tmpdecor];
        end 
        decor = [decorj 99 decor];
        objr_whole = [objr 99 objr_whole];
        coef_whole = [coef 99*ones([ncov 1]) coef_whole];
        lambdas_whole = [lambdas 99 lambdas_whole];

    end
end


