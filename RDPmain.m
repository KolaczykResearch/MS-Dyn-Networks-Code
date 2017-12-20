%% Simulation study using the recursive dyadic partition
%  There are three models (Model A, Model B, Model C) in the paper. This
%  part is for Model B. For model A and C, the only difference is how the
%  data is simulated.
%% Result:  count    := equals 1 if change point detected
%           coeflist := list of vectors of coefficients, each level is
%           separated by a vector of 99.
%           decorlist:= list of vectors of split points.
%% Calls:   Requires the cvx package
%  More information can be found in the paper linked at:
%  http://www.stanford.edu/~boyd/papers/distr_opt_stat_learning_admm.html
cvx_setup  
randn('seed', 0);
rand('seed',0);
coeflist = {};
decorlist= {};
count = [];

iter = 100;
K = 3;
p = 2;
parpool(4);

% Two change points; one in 1/2, the other one in 3/4
parfor i = 1:iter,
T = 514;

A = zeros(K, K, p);
A(1, 2, 1) = 0.5;
A(1, 2, 2) = 0.25;
Data1 = data_gen(T, A);

T = 257;
A = zeros(K, K, p);
A(1, 3, 1) = 0.5;
A(1, 3, 2) = 0.25;
Data2 = data_gen(T, A);

A = zeros(K, K, p);
A(1, 2, 1) = -0.5;
A(1, 3, 1) = -0.5;
Data3 = data_gen(T, A);

Data = [Data1; Data2; Data3];

y_id = 1;
[bestFit, decor, lambda, coef, objr] = rdp_ts(y_id, Data, p, 0.05);
count(i) = decor(1);
count2(i) = decor(4);
decorlist{i} = decor;
coeflist{i} = coef;
end

