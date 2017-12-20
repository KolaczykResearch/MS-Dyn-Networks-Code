%% Simulation study using the recursive partition
%  There are three models (Model A, Model B, Model C) in the paper. This
%  part is for Model B. For model A and C, the only difference is how the
%  data is simulated.
%% Result:  bplist := list of change points
%% Calls:   Requires the cvx package
%  More information can be found in the paper linked at:
%  http://www.stanford.edu/~boyd/papers/distr_opt_stat_learning_admm.html

cvx_setup
randn('seed', 0);
rand('seed',0);
bplist = {};


K = 3;
p = 2;
iter = 100;
parpool(4);

parfor i = 1:iter
T = 512;
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
[Y, X] = covar_design(Data, p, y_id);

block = repmat(p, K, 1);
block(y_id) = 0;

nfeature = size(X, 2);
len = size(X, 1);
bp = [];

[t0, t1, mm] = TriangleMDN(y_id, Data, p, 0.05);
bp = GetChangepoints(bp, 1, len, t0, t1, mm);
bp = sort(bp);

bplist{i} = bp;
end

