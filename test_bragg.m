%% Example bragg grating 
% create a Brgg grating that provides 30dB of suppresion at 980nm
% use materials with ref index of n1 = 1.5 and n2 = 2
% 
% for ref index of 1.5 and 2 
% d = lambda0/4n therefor 
% d1 = 980nm/(4*1.5) = 163nm
% d2 = 980nm/(4*2.0) = 122nm
% assume ur = 1
% n = sqrt(e_r*u_r)
% References: 
% https://en.wikipedia.org/wiki/Fresnel_equations
% https://en.wikipedia.org/wiki/Refractive_index

close all

addpath('./src')

param = getDefaultParameters();
param.nzsrc = 10;
param.spacerRegion = 50;
param.num_plot = 10;

lambda = 800e-9:0.1e-9:1200e-9;
f_Hz = flip(param.const.c0./lambda);

n1 = 1.5;
n2 = 2;
d1 = 980e-9/(4*n1);
d2 = 980e-9/(4*n2);

lambda0 = 980e-9;
dmin = min(d1,d2);

param = get_resolution(param, lambda, dmin, max(n1,n2));
% calculate total number of points for simulation
numRep = 40;
param.Nz = numRep*round(d2./param.dz)  + numRep*round(d1./param.dz) + 2*param.spacerRegion + 3; % 3 comes from T R and source cells

% material properties
ER = ones(1,param.Nz);
UR = ones(1,param.Nz);

% calculate start and end points of slab
nz_end = param.spacerRegion + 2; % end of spacer region, start of grating
d1numCells = round(d1./param.dz);
d2numCells = round(d2./param.dz);

for n = 1:numRep
    nz_start = nz_end;
    nz_end = nz_start+d2numCells;
    ER(nz_start:nz_end) = n2^2;

    nz_start = nz_end;
    nz_end = nz_start+d1numCells;
    ER(nz_start:nz_end) = n1^2;
end

% create ref index structure
n_matrix = sqrt(UR.*ER);

figure; stem(n_matrix);

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