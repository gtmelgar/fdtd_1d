lambda = 200:300;
n_max = 1;
lambda_res = 10;
c0 = 299792458; % free space light speed m/s

dmin = 100;
nd_res = 4;
n_bc = 1;

% get default resolution
dz1 = min(lambda)/n_max/lambda_res;
dz2 = dmin/nd_res;
dz = min(dz1,dz2);

dt = n_bc.*dz/(2*c0);

dc = 100;

% further refine dz
N = ceil(dc/dz);
dz = dc/N;

Nz = 500;

% get num of time steps to solve for
t_prop = n_max*Nz*dz./c0;
pulse_duration = 100;

num_time_steps = 15.*pulse_duration + 5 * t_prop;

% generate gaussian pulse and delay by 6tau
fmax = c0/min(lambda);
tau = 0.5./fmax;
t = 0:dt:num_time_steps*dt;
g_t = exp(-(t-6*tau).^2./tau.^2); 

% init FFT 
num_steps = 1001;
freq = linspace(-fmax/2,fmax/2,num_steps);

num_f = num_steps;
K = exp(-1i*2*pi*dt*freq);
EyR = zeros(1,num_f);
EyT = zeros(1,num_f);
src = zeros(1,num_f);

% material properties
ER = ones(1,Nz);
UR = ones(1,Nz);

% update coeffs
mEy = (c0*dt)./ER;
mHx = (c0*dt)./UR;

Hx = zeros(1,Nz);
Ey = zeros(1,Nz);

E2 = 0;
E1 = 0;
H2 = 0;
H1 = 0;

%location of source
nzsrc = 100;

spacer_region = max(lambda);
% solve
for T = 1:num_time_steps
    
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

    Ey(nzsrc) = Ey(nzsrc) + g_t(T);

    % update fourier transforms
    % Update Fourier Transforms
    % for nf = 1 : num_f
    %     EyR(nf) = EyR(nf) + (K(nf)^T)*Ey(1);
    %     EyT(nf) = EyT(nf) + (K(nf)^T)*Ey(Nz);
    %     src(nf) = src(nf) + (K(nf)^T)*g_t(Nz);
    % end
    if mod(T,2000)
        clf;
        plot(Ey);hold on 
        plot(Hx);
        ylim([-1.5, 1.5])
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
