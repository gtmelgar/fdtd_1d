function [Esrc,Hsrc,UR,ER] = gen_pulse(ER, UR, n_matrix, param, tau, time_vector, lambda)

t0 = 6*tau;

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