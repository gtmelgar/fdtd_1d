% example simulaition 
% 1foot slab (30.48cm) ur = 1 , e_r = 12 with air around it
% apply AR coating to allow 2.4GHz transmission

% Assume no magnetic response 
% AR coating thickness should be 1/4 wavelength

close all
clear 
addpath('./src')
param = getDefaultParameters();
param.spacerRegion = 100;
param.num_plot = 100;
param.num_bins = 500;
f_Hz = (0.1:0.01:5)*1e9;
lambda = param.const.c0./f_Hz;

er_1 = 12;
er_air = 1;

n1 = sqrt(er_1*1);

% calculate geometric mean
er_2 = sqrt(er_1*er_air); 

% calculate ref. index
n2 = sqrt(er_2);

lambda0 = param.const.c0./2.4e9;
d2 = lambda0./(4*n2);


dmin = min(30.48e-2,d2); % slab to simulate
[param] = get_resolution(param, lambda, dmin, max(n1,n2));
% calculate total number of points for simulation
param.Nz = param.N*2 + round(30.48e-2./param.dz) + 2*param.spacerRegion + 3; % 3 comes from T R and source cells

% material properties
ER = ones(1,param.Nz);
UR = ones(1,param.Nz);

% calculate start and end points of slab
Nz1 = param.spacerRegion + 2; % start of AR coating
Nz2 = Nz1 + param.N; % end of AR coating
Nz3 = Nz2 + round(30.48e-2./param.dz); % material + start 2nd set of AR coating
Nz4 = Nz3 + param.N; %end 2nd set of AR coating

% create structure
UR(Nz1:Nz2) = 1;
ER(Nz1:Nz2) = er_2;

UR(Nz2:Nz3) = 1;
ER(Nz2:Nz3) = er_1;

UR(Nz3:Nz4) = 1;
ER(Nz3:Nz4) = er_2;

% create ref index structure
n_matrix = sqrt(UR.*ER);

% compute min time step
param.dt = min(n_matrix).*param.dz/(2*param.const.c0);


% generate gaussian pulse and delay by 6tau
fmax = param.const.c0/min(lambda);
tau = 0.5./fmax;

% get num of time steps to solve for
t_prop = max(n_matrix)*param.Nz*param.dz./param.const.c0;

sim_time = 12.*tau + 5 * t_prop;
num_steps = ceil(sim_time./param.dt);

time_vector = 0:param.dt:(num_steps-1)*param.dt;

[Esrc,Hsrc,UR,ER] = gen_pulse(ER, UR, n_matrix, param, tau, time_vector, lambda);

run_fdtd(ER,UR,Esrc,Hsrc,num_steps,param,f_Hz)
