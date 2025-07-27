function run_fdtd(ER,UR,Esrc,Hsrc,num_steps,param,f_Hz)

n_matrix = sqrt(UR.*ER);

% init FFT 
EyR = zeros(1,param.num_bins);
EyT = zeros(1,param.num_bins);
src = zeros(1,param.num_bins);

f_max_plot = max(f_Hz);
freq = linspace(0,f_max_plot,param.num_bins);
K = exp(-1i*2*pi*param.dt*freq);


% update coeffs
mEy = (param.const.c0*param.dt)./ER;
mHx = (param.const.c0*param.dt)./UR;

Hx = zeros(1,param.Nz);
Ey = zeros(1,param.Nz);

E2 = 0;
E1 = 0;
H2 = 0;
H1 = 0;

param.do_plot = 1;

% solve
for T = 1:num_steps
    
    % H2 and H1 are used for perfect boundary condition
    % H2 and H1 are needed since wave propragates 1 cell every two time steps
    H2=H1; 
    H1=Hx(1);

    % update H from E
    for nz = 1:param.Nz-1
        Hx(nz) = Hx(nz) + mHx(nz)*(Ey(nz+1) - Ey(nz))/param.dz;
    end
    % note E2 = 0 are Dirichlet boundary conditions 
    % for periodic boundary conditions use E2 = Ey(1) 
    Hx(param.Nz) = Hx(param.Nz) + mHx(param.Nz)*(E2 - Ey(param.Nz))/param.dz; 
    
    % TF/SF correction term
    Hx(param.nzsrc-1) = Hx(param.nzsrc-1) - mHx(param.nzsrc-1) .* Esrc(T) ./ param.dz;

    % E2 and E1 are used for perfect boundary condition
    % E2 and E1 are needed since wave propragates 1 cell every two time steps
    E2=E1; 
    E1=Ey(param.Nz);

    % note H2 = 0 are Dirichlet boundary conditions
    % for periodic boundary conditions use H2 = Hx(Nz) 
    Ey(1) = Ey(1) + mEy(1)*(Hx(1) - H2)/param.dz;

    % update E from H 
    for nz = 2:param.Nz
        Ey(nz) = Ey(nz) + mEy(nz)*(Hx(nz) - Hx(nz-1))/param.dz;
    end
    % Ey(param.nzsrc) = Ey(param.nzsrc) + g_t(T); % soft source for debug only

    % TF/SF correction term
    Ey(param.nzsrc) = Ey(param.nzsrc) - mEy(param.nzsrc) * Hsrc(T) ./ param.dz;
    
    % update fourier transforms
    % Update Fourier Transforms
    for nf = 1 : param.num_bins
        EyR(nf) = EyR(nf) + (K(nf)^T)*Ey(1);
        EyT(nf) = EyT(nf) + (K(nf)^T)*Ey(param.Nz);
        src(nf) = src(nf) + (K(nf)^T)*Esrc(T);
    end

    if mod(T,param.num_plot)==0 && param.do_plot
        
        clf;
        subplot(3,1,1)
        plot(Ey);hold on 
        plot(Hx);
        ylim([-1.5, 1.5])
        yyaxis right
        plot(n_matrix);
        grid on;
        
        subplot(3,1,2)
        plot(abs(EyR./src).^2);hold on
        plot(abs(EyT./src).^2);
        plot(abs(EyT./src).^2+abs(EyR./src).^2);
        ylim([-0.5, 1.5])

        subplot(3,1,3)
        plot(abs(EyR)); hold on
        plot(abs(EyT))
        plot(abs(src))
        drawnow;
    end
end

% scale FFT
EyR = EyR*param.dt;
EyT = EyT*param.dt;

% transmit and reflectance plots
reflectance = abs(EyR./src).^2;
transmittance = abs(EyT./src).^2;
CON = reflectance + transmittance;
