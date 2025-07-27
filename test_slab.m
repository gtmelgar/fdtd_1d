% example simulaition 
% 1foot slab (30.48cm) ur = 2 , e_r = 6 with air around it
% test R and T from 0 to 1GHz
close all
tic;

f_Hz = (0.1:0.1:1)*1e9;

u_r = 2;
e_r = 6;

n_max = sqrt(u_r*e_r);

lambda_res = 20;
c0 = 299792458; % free space light speed m/s

lambda = c0./f_Hz;

dmin = 30.48e-2;
nd_res = 4;

% get default resolution
dz1 = min(lambda)/n_max/lambda_res;
dz2 = dmin/nd_res;
dz = min(dz1,dz2);

dc = dmin;

% further refine dz (N is num cells to fit device)
N = ceil(dc/dz);
dz = dc/N; 

% calculate num of grid points
% spacerRegion = (ceil(max(lambda)/dz)); %minumum for 2D/3D
spacerRegion = 1000; % spacer not neede for 1D
Nz = N + 2*spacerRegion + 3; % 3 comes from T R and source cells

Nz1 = spacerRegion + 2;
Nz2 = Nz1 + N - 1;

% material properties
ER = ones(1,Nz);
UR = ones(1,Nz);

UR(Nz1:Nz2) = u_r;
ER(Nz1:Nz2) = e_r;

n_matrix = sqrt(UR.*ER);

% compute min time step
dt = min(n_matrix).*dz/(2*c0);

% generate gaussian pulse and delay by 6tau
fmax = c0/min(lambda);
tau = 0.5./fmax;
t0 = 6*tau;

% get num of time steps to solve for
t_prop = n_max*Nz*dz./c0;

sim_time = 12.*tau + 5 * t_prop;
num_steps = ceil(sim_time./dt);

time_vector = 0:dt:(num_steps-1)*dt;

% TF/SF formulation
nzsrc = 2;
ersrc = ER(nzsrc); % TODO: get from the source injection point
ursrc = UR(nzsrc);
nsrc = n_matrix(nzsrc);

A = -sqrt(ersrc/ursrc);
% delay between E and H
delt = nsrc*dz/(2*c0) + dt/2;

g_t = exp(-((time_vector-t0)./tau).^2); % e-field
Esrc = g_t;
Hsrc = A*exp(-((time_vector-t0+delt)/tau).^2);

% init FFT 
num_bins = 101;
f_max_plot = 1e9;
freq = linspace(0,f_max_plot,num_bins);
K = exp(-1i*2*pi*dt*freq);
EyR = zeros(1,num_bins);
EyT = zeros(1,num_bins);
src = zeros(1,num_bins);

% numerical disersion correction factor
n_avg = mean(n_matrix);
f_mid = mean(c0./lambda);
k0 = 2*pi*f_mid./c0;
disp_f = c0*dt/(n_avg*dz) * sin(k0*n_avg*dz/2)/sin(c0*k0*dt/2);
UR = disp_f*UR;
ER = disp_f*ER;

% update coeffs
mEy = (c0*dt)./ER;
mHx = (c0*dt)./UR;

Hx = zeros(1,Nz);
Ey = zeros(1,Nz);

E2 = 0;
E1 = 0;
H2 = 0;
H1 = 0;

% %location of source
% nzsrc = 100;
do_plot = 1;

% solve
for T = 1:num_steps
    
    % H2 and H1 are used for perfect boundary condition
    % H2 and H1 are needed since wave propragates 1 cell every two time steps
    H2=H1; 
    H1=Hx(1);

    % update H from E
    for nz = 1:Nz-1
        Hx(nz) = Hx(nz) + mHx(nz)*(Ey(nz+1) - Ey(nz))/dz;
    end
    % note E2 = 0 are Dirichlet boundary conditions 
    % for periodic boundary conditions use E2 = Ey(1) 
    Hx(Nz) = Hx(Nz) + mHx(Nz)*(E2 - Ey(Nz))/dz; 
    
    % TF/SF correction term
    Hx(nzsrc-1) = Hx(nzsrc-1) - mHx(nzsrc-1) .* Esrc(T) ./ dz;

    % E2 and E1 are used for perfect boundary condition
    % E2 and E1 are needed since wave propragates 1 cell every two time steps
    E2=E1; 
    E1=Ey(Nz);

    % note H2 = 0 are Dirichlet boundary conditions
    % for periodic boundary conditions use H2 = Hx(Nz) 
    Ey(1) = Ey(1) + mEy(1)*(Hx(1) - H2)/dz;

    % update E from H 
    for nz = 2:Nz
        Ey(nz) = Ey(nz) + mEy(nz)*(Hx(nz) - Hx(nz-1))/dz;
    end
    % Ey(nzsrc) = Ey(nzsrc) + g_t(T); % soft source for debug only

    % TF/SF correction term
    Ey(nzsrc) = Ey(nzsrc) - mEy(nzsrc) * Hsrc(T) ./ dz;
    
    % update fourier transforms
    % Update Fourier Transforms
    for nf = 1 : num_bins
        EyR(nf) = EyR(nf) + (K(nf)^T)*Ey(1);
        EyT(nf) = EyT(nf) + (K(nf)^T)*Ey(Nz);
        src(nf) = src(nf) + (K(nf)^T)*Esrc(T);
    end

    if mod(T,10)==0 && do_plot
        
        clf;
        subplot(2,1,1)
        plot(Ey);hold on 
        plot(Hx);
        plot(UR/max(UR));
        grid on;
        ylim([-1.5, 1.5])
        subplot(2,1,2)
        plot(abs(EyR./src).^2);hold on
        plot(abs(EyT./src).^2);
        plot(abs(EyT./src).^2+abs(EyR./src).^2);
        ylim([-0.5, 1.5])
        drawnow;
    end
end

% scale FFT
EyR = EyR*dt;
EyT = EyT*dt;

% transmit and reflectance plots
reflectance = abs(EyR./src).^2;
transmittance = abs(EyT./src).^2;
CON = reflectance + transmittance;
toc