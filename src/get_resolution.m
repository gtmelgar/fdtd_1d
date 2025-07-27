function [param] = get_resolution(param, lambda, dmin, n_max)

% get default resolution
dz1 = min(lambda)/n_max/param.lambda_res;
dz2 = dmin/param.nd_res;
param.dz = min(dz1,dz2);

dc = dmin;

% further refine dz (N is num cells to fit device)
param.N = ceil(dc/param.dz);
param.dz = dc/param.N; 

% calculate total number of points for simulation
param.Nz = param.N + 2*param.spacerRegion + 3; % 3 comes from T R and source cells

