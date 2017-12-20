%% MEG analysis
%  We first import the break points detected using MEG_breakpoints.m.
%% Result:  This file generate all plots of the application part of the 
%  ?Dynamic Networks with Multi-scale Temporal Structure? paper?
%% Calls:   Requires the cvx package
%  More information can be found in the paper linked at:
%  http://www.stanford.edu/~boyd/papers/distr_opt_stat_learning_admm.html
% Print the change point of the two areas 
load bp_fp.mat
load bp_vs.mat

plotbp(store_bp_vs,'black');
xlabel('Time(ms)');
ylabel('Count');

plotbp(store_bp_fp, 'black');
xlabel('Time(ms)');
ylabel('Count');

% VS network
coef_vs = {};
for i = dt,
    Data = out_dat(:, :, i).';
    ds = mat2dataset(Data, 'VarNames', {'IPS_R','MPFC_R','MTp_R','SPL_R','FEF_R','STP_R','VIP_R','V3a_R','aud_R','IPS_L','MPFC_L','MTp_L','SPL_L','FEF_L','STP_L','VIP_L','V3a_L','aud_L'});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Extract ROI Information %%%%%%%%%%%%%%%%%%%% 
    rh_circle = double(ds(:, [8 3 7 4 5 2])); % Two networks: V3a, MT+, VIP; and SPL, FEF, MPFC
    
    current_split = 0;
    Data = rh_circle;
    
    K = size(Data, 2);
    bp = sort(unique(store_bp_vs{i}));
    
    coef_vs{i} = RecoverCoef(Data, bp, 0.05, p);
end

test_coef = cell2mat(coef_vs);
coef_vs_mat = reshape(test_coef, [6, 6, 160, 1499]); % row, column, trial, time

shadedplot(coef_vs_mat(1,2,:,:), 'MT+     --->     V3a')
shadedplot(coef_vs_mat(2,1,:,:), 'V3a     --->     MT+')
shadedplot(coef_vs_mat(1,3,:,:), 'V3a     --->     VIP')
shadedplot(coef_vs_mat(3,1,:,:), 'VIP     --->     V3a')
shadedplot(coef_vs_mat(2,3,:,:), 'MT+     --->     VIP')
shadedplot(coef_vs_mat(3,2,:,:), 'VIP     --->     MT+')

% FP network
% VS network
coef_fp = {};
for i = dt,
    Data = out_dat(:, :, i).';
    ds = mat2dataset(Data, 'VarNames', {'IPS_R','MPFC_R','MTp_R','SPL_R','FEF_R','STP_R','VIP_R','V3a_R','aud_R','IPS_L','MPFC_L','MTp_L','SPL_L','FEF_L','STP_L','VIP_L','V3a_L','aud_L'});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Extract ROI Information %%%%%%%%%%%%%%%%%%%% 
    rh_circle = double(ds(:, [8 3 7 4 5 2])); % Two networks: V3a, MT+, VIP; and SPL, FEF, MPFC
    
    current_split = 0;
    Data = rh_circle;
    
    K = size(Data, 2);
    bp = sort(unique(store_bp_fp{i}));
    
    coef_fp{i} = RecoverCoef(Data, bp, 0.05, p);
end

test_coef = cell2mat(coef_fp);
coef_fp_mat = reshape(test_coef, [6, 6, 160, 1499]); % row, column, trial, time

shadedplot(coef_fp_mat(4,5,:,:), 'SPL     --->     FEF')
shadedplot(coef_fp_mat(5,4,:,:), 'FEF     --->     SPL')
shadedplot(coef_fp_mat(4,6,:,:), 'SPL     --->     DLPFC')
shadedplot(coef_fp_mat(6,4,:,:), 'DLPFC     --->     SPL')
shadedplot(coef_fp_mat(5,6,:,:), 'FEF     --->     DLPFC')
shadedplot(coef_fp_mat(6,5,:,:), 'DLPFC     --->     FEF')
