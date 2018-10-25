%% Numerical Solution to Quasi-1D Non-Conservation form Nozzle flow 
%% using MacCormak's technique

clear all
close all
clc


% Inputs
% n = input('No of grid points: ');
n = 31;
L = 3;
x = linspace(0,L,n);
dx = L/(n-1);
gamma = 1.4;


% Initial profiles

rho = 1-0.3146*x;
t = 1-0.2314*x; % t = temp
v = (0.1+1.09*x).*t.^0.5;

% Area
a = 1+2.2*(x - 1.5).^2;

% Time steps
nt = 1400;
% C = input('Input the required Courant no: ');
C = 0.5;
dt = min((C.*dx)./(sqrt(t)+v));

% Outer time loop
tic
for k = 1:nt
    
    rho_old = rho;
    v_old = v;
    t_old = t;
    
    
    %% Predictor method
   
    for j = 2:n-1
        
        % Predictor method Governing equations
        
        [drho_dt_p(j),dvdt_p(j),dtdt_p(j)] = non_conservation_p(v,rho,t,a,j,dx,gamma);
        
        %Solution update
        rho(j) = rho(j) + drho_dt_p(j)*dt;
        v(j) = v(j) + dvdt_p(j)*dt;
        t(j) = t(j) + dtdt_p(j)*dt;
    end
    
    %% Corrector method
    for j = 2:n-1
        
        % Corrector method Governing equations
        
        [drho_dt_c(j),dvdt_c(j),dtdt_c(j)] = non_conservation_c(v,rho,t,a,j,dx,gamma);
        
    end
    
    % Compute the average time derivative
    
    dvdt = 0.5*(dvdt_p + dvdt_c);
    drhodt = 0.5*(drho_dt_p + drho_dt_c);
    dtdt = 0.5*(dtdt_p+dtdt_c);
    
    % Solution update
    for i = 2:n-1
        v(i) = v_old(i) + dvdt(i)*dt;
        rho(i) = rho_old(i) + drhodt(i)*dt;
        t(i) = t_old(i) + dtdt(i)*dt;
    end
    
    % Boundary conditions
    
    % Inlet B.C.
    v(1) = 2*v(2) - v(3);
    
    %Outlet B.C.
    v(n) = 2*v(n-1) - v(n-2);
    rho(n) = 2*rho(n-1) - rho(n-2);
    t(n) = 2*t(n-1) - t(n-2);
    
    m = rho.*v.*a;
    mach = v./(sqrt(t));
    P = rho.*t;
    
    % Solution at throat (area = 1)
    rho_th(k) = rho(16);
    temp_th(k) = t(16);
    Pr_th(k) = P(16);
    mach_th(k) = mach(16);
    v_th(k) = v(16);
    m_th(k) = m(16);
    k_th = 1:nt;
end

% figure(1)
% plot(k_th,m_th,'-m'); hold on;
% legend('Non-conservation');

figure

subplot(2,2,1);
plot(x,v); grid on;
xlabel('x/L');
ylabel('v/a_o');
title('Velocity profile');

subplot(2,2,2);
plot(x,mach); grid on;
xlabel('x/L');
ylabel('Mach no m');
title('Mach no profile');

subplot(2,2,3);
plot(x,rho); grid on;
xlabel('x/L');
ylabel('{\rho_e}/{\rho_o}');
title('Density profile');

subplot(2,2,4);
plot(x,m); grid on;
xlabel('x/L');
ylabel('{\rho}av/{\rho_o}a_ov_o');
title('Mass flow rate profile');

suptitle('Parameter profile for 1D nozzle flow for Non-Conservation form');

%% Writing output in word file
fid = fopen('nonconservation.doc ','w');
fprintf(fid,'Parameter for 1D Conservation form nozzle flow \n');
fprintf(fid,'for grid size %d and Courant no %f \n\n', n,C);
fprintf(fid,'Velocity  Density  Temp   Mach no  Massflowrate \n\n');
fprintf(fid, '%f %f %f %f %f \n', [v;rho;t;mach;m]);
fclose(fid);

%plot(x,m);
toc
    
    
    