% example simulaition 
% 1foot slab (30.48cm) ur = 2 , e_r = 6 with air around it
% test R and T from 0 to 1GHz
close all
tic;

addpath('./src')

param = getDefaultParameters();
param.nzsrc = 100;
param.spacerRegion = 200;
f_Hz = (0.1:0.1:1)*1e9;

u_r = 2;
e_r = 6;
n_max = sqrt(u_r*e_r);

lambda = param.const.c0./f_Hz;

dmin = 30.48e-2; % slab to simulate

[param] = get_resolution(param, lambda, dmin, n_max);

% material properties
ER = ones(1,param.Nz);
UR = ones(1,param.Nz);

% calculate start and end points of slab
Nz1 = param.spacerRegion + 2;
Nz2 = Nz1 + param.N - 1;

% create structure
UR(Nz1:Nz2) = u_r;
ER(Nz1:Nz2) = e_r;

% create ref index structure
n_matrix = sqrt(UR.*ER);

% compute min time step
param.dt = min(n_matrix).*param.dz/(2*param.const.c0);

% generate gaussian pulse and delay by 6tau
fmax = param.const.c0/min(lambda);
tau = 0.5./fmax;
t0 = 6*tau;

% get num of time steps to solve for
t_prop = n_max*param.Nz*param.dz./param.const.c0;

sim_time = 12.*tau + 5 * t_prop;
num_steps = ceil(sim_time./param.dt);

time_vector = 0:param.dt:(num_steps-1)*param.dt;

% TF/SF formulation
ersrc = ER(param.nzsrc);
ursrc = UR(param.nzsrc);
nsrc = n_matrix(param.nzsrc);

A = -sqrt(ersrc/ursrc);
% delay between E and H
delt = nsrc*param.dz/(2*param.const.c0) + param.dt/2;
g_t = exp(-((time_vector-t0)./tau).^2); % e-field
Esrc = g_t;
Hsrc = A*exp(-((time_vector-t0+delt)/tau).^2);

% numerical disersion correction factor
n_avg = mean(n_matrix);
f_mid = mean(param.const.c0./lambda);
k0 = 2*pi*f_mid./param.const.c0;
disp_f = param.const.c0*param.dt/(n_avg*param.dz) * sin(k0*n_avg*param.dz/2)/sin(param.const.c0*k0*param.dt/2);
UR = disp_f*UR;
ER = disp_f*ER;

run_fdtd(ER,UR,Esrc,Hsrc,num_steps,param,f_Hz)

toc