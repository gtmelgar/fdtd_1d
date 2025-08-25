function run_fdtd(ER,UR,Esrc,Hsrc,num_steps,param,f_Hz)

n_matrix = sqrt(UR.*ER);
num_bins = length(f_Hz);
% init FFT 
EyR = zeros(1,num_bins);
EyT = zeros(1,num_bins);
src = zeros(1,num_bins);

Esrc_FFT = abs(fft(Esrc));
Esrc_FFT_x = (0:length(Esrc)-1)*1./param.dt*1/length(Esrc);
% f_max_plot = max(f_Hz);
% freq = linspace(0,f_max_plot,num_bins);
freq = f_Hz;
K = exp(-1i*2*pi*param.dt*f_Hz);


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
    for nf = 1 : length(f_Hz)%num_bins
        EyR(nf) = EyR(nf) + (K(nf)^T)*Ey(1);
        EyT(nf) = EyT(nf) + (K(nf)^T)*Ey(param.Nz);
        src(nf) = src(nf) + (K(nf)^T)*Esrc(T);
    end

    if mod(T,param.num_plot)==0 && param.do_plot
        
        clf;
        subplot(3,1,1)
        plot(Ey);hold on 
        plot(Hx);
        xlabel('unit cell')
        ylabel('amplitude')
        ylim([-1.5, 1.5])
        yyaxis right
        plot(n_matrix);
        ylabel('refractive index')
        % legend({'Ey','Hx','Ref Index'})
        % grid on;
        title(['Field at step ' num2str(T) ' of ' num2str(num_steps)])
        
        subplot(3,1,2)
        plot(freq,abs(EyT./src).^2);hold on
        plot(freq,abs(EyR./src).^2);
        plot(freq,abs(EyT./src).^2+abs(EyR./src).^2);
        ylim([-0.5, 1.5])

        % semilogy(freq,abs(EyT./src).^2);hold on
        % semilogy(freq,abs(EyR./src).^2);
        % semilogy(freq,abs(EyT./src).^2+abs(EyR./src).^2);

        % plot(freq,10*log10(abs(EyT./src).^2));hold on
        % plot(freq,10*log10(abs(EyR./src).^2));
        % plot(freq,10*log10(abs(EyT./src).^2+abs(EyR./src).^2));
        % ylim([-30 1])

        
        % set(gca, 'YScale', 'log') 
        % % legend({'Transmittance','Reflectance','T + R'})
        title(['Field at step ' num2str(T*param.dt/1e-9) 'ns of ' num2str(num_steps*param.dt/1e-9) 'ns'])

        subplot(3,1,3)
        plot(Esrc_FFT_x,Esrc_FFT);hold on
        plot(freq,abs(EyT));
        plot(freq,abs(EyR))
        plot(freq,abs(src))
        xlim([freq(1) freq(end)]);
        % legend({'T Spectrum','R Spectrum','Source Spectrum'})

        % set(gca,'LegendColorbarListeners',[]);
        % setappdata(gca,'LegendColorbarManualSpace',1);
        % setappdata(gca,'LegendColorbarReclaimSpace',1);

        drawnow;
    end
end

% % % scale FFT
% % EyR = EyR*param.dt;
% % EyT = EyT*param.dt;
% % 
% % % transmit and reflectance plots
% % reflectance = abs(EyR./src).^2;
% % transmittance = abs(EyT./src).^2;
% % CON = reflectance + transmittance;
