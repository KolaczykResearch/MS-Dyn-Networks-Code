%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Name:    covar_design.m        Created: 11/14/16    Revised: 
%
%% Usage:   Reshape the data into a form that suitable for group-lasso
%           regression. The covariate matrix looks like:
%           data looks like [x1,2 x1,3 x2,2 x2,3...]
%                           [x1,3 x1,4 x2,3 x2,4...]
%                           [x1,t x1,t+1 ..........]
%           Removes the oldest p obs from the response variable and removes
%           the newest p obs from the independent variables.
%
%% Inputs:  X    := Time series data
%           p    := Number of lags used in the VAR model
%           y_id := Index of the series which will be used as the response
%                   variable
%
%%  Output:  Y            := Response variable
%           Covariates    := Design matrix 
%% Calls:   Only internal Matlab functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Y, covariates] = covar_design(X, p, y_id)

[T, n] = size(X); % size of the network
mat = [];
ind = [];

% We store the lag p series in the pth layer of the 3-dimension array mat
% Return a matrix with the first i rows removed (The oldest is getting removed)
% This array will be reshaped later as an appropriate design matrix
for i = 1:p
    mat(:, :, i) = lagmatrix(X, i); 
end

% first b columns are all independent variables of lag 1, next b columns are of lag 2...
[a, b, c] = size(mat);
tmp = reshape(mat, [a b*c]);

for i = 1:n
    ind = [ind i:n:i+(p-1)*n];
end

ind = ind.';
data = tmp(:, ind);

covariates = removerows(data, 'ind', 1:p);
Y = X((1+p):T, y_id); 
covariates(:, (p*y_id-p+1):(p*y_id)) = [];  % keep this if self loop not allowed

end