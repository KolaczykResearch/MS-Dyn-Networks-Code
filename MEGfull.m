%% Brain image analysis
%  We extract the regions of interest and use the partition based
%  multiscale causal network work model to detect the changepoints among
%  the regions.
%% Result:  store_bp: break(change) points detected from the network (regions of interest)
%% Calls:   Requires the cvx package
%  More information can be found in the paper linked at:
%  http://www.stanford.edu/~boyd/papers/distr_opt_stat_learning_admm.html
% tic
cvx_setup
tic
path = '/Users/Xinyu/Documents/BU/Research/Multiscale Model/Data/alpha7/datafile/';
d = dir(path);
names = {d.name};

f = 3;  %
filename = strcat(path, names(f));
importfile(char(filename));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = size(out_dat, 3);                  % number of trials for the object
p = 7;                                 % lags used in the granger causal model
store_bp = {};
alpha = 0.05;                          % control the probability of making a type I error

parpool(8);
parfor i = 1:n
    Data = out_dat(:, :, i).';
    ds = mat2dataset(Data, 'VarNames', {'IPS_R','DLPFC_R','MTp_R','SPL_R','FEF_R','STP_R','VIP_R','V3a_R','aud_R','IPS_L','DLPFC_L','MTp_L','SPL_L','FEF_L','STP_L','VIP_L','V3a_L','aud_L'});
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Extract ROI Information %%%%%%%%%%%%%%%%%%%%

     rh_circle = double(ds(:, [8 3 7 4 5 2])); % Two networks: V3a, MT+, VIP; and SPL, FEF, DLPFC
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    data = diff(rh_circle);
    K = size(rh_circle, 2); % size of the network
    current_split = 0;
    
%    data = removeAR(data, AR, MA);
     data = double(rh_circle);

    for y_id = 1:6 % loop through each node in the network
        len = size(data, 1) - p;
        bp = [];

        [t0, t1, mm] = TriangleMDN(y_id, double(data), p, alpha);
        bp = GetChangepoints(bp, 1, len, t0, t1, mm);
        store_bp{i} = bp;
    end
end
save store_bp 
toc

