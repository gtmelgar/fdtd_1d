function [param] = get_resolution(param, lambda, dmin, n_max)

% get default resolution
dz1 = min(lambda)/n_max/param.lambda_res;
dz2 = dmin/param.nd_res;
param.dz = min(dz1,dz2);

dc = dmin;

% further refine dz (N is num cells to fit smallest device length)
param.N = ceil(dc/param.dz);
param.dz = dc/param.N; 

