function param = getDefaultParameters()

param.const.c0 = 299792458; % speed of light m/s

param.lambda_res    = 20;  % number of points to resolve one wavelength
param.nd_res        = 4;   % number of points to resolve smallest feature
param.spacerRegion  = 500; % not needed for 1D, at least 1lambda for 2D/3D

param.nzsrc         = 2;   % source injection point
param.num_bins      = 100; % num bins for for FFT

param.num_plot      = 10;  % plot every n time steps
param.do_plot       = 1;   % enable (1), disable (0) plots 

% grid properties
param.dz    = 0; % smallest grid resolution
param.N     = 0; % number of cells in structure
param.Nz    = 0; % number of total cell in simulation

% time properties
param.dt    = 0; % time step resolution