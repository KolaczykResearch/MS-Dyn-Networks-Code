%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Name:    data_gen.m        Created: 11/06/16    Revised: 11/14/16
%
%% Usage:   Generate multi-variate time series data of length T using
%           coeffcient matrix A
%
%% Inputs:  T := length of the time series
%           A := an [K K p] dimensional coeffcient matrix A where K is the
%           number of time series and p is the number of lags in the VAR
%           model
%
%% Output:  X := T by K dimensional vector autoregressive process.
%% Calls:   Only internal Matlab functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X] = data_gen(T, A)

[K, ~, p] = size(A);
epsilon = eye(K);
C = cell(1, p);

for i = 1:p
    C{i} = A(:,:,i);
end

struct = vgxset('n', K, 'AR', C, 'Q', epsilon);
X = vgxsim(struct, T);

