function [data_scaled] = fn_scaling(data,mean_target,max_target)
% Exponential scaling
% Should not allow negative data
% loadfactor = mean/max

hf_datascale = @(data,K,F) K.*data.^F;

optfun = @(X) [...
    mean_target - mean(hf_datascale(data,X(1),X(2)))...
    max_target - max(hf_datascale(data,X(1),X(2)))...
    ];

X = fsolve(optfun,[1 1]);

data_scaled = hf_datascale(data,X(1),X(2));