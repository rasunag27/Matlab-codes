% Solving to Quasi 1D nozzle flow using MacCormak technique
% with the usage of Conservation form of Governing equation

clear all
close all
clc

% Initialization
% n = input('No of grid points: ');
n = 31;
L = 3;
x = linspace(0,L,n);
dx = L/(n-1);
gamma = 1.4;

%% Initial profiles and Conditions

% Area

a = 1+2.2*(x - 1.5).^2;

%Density and temp Initial conditions

rho = zeros(1,n);
t = zeros(1,n);
for i = 1:n
    if (x(i)>=0 && x(i)<=0.5)
        rho(i) = 1;
        t(i) = 1;
    elseif (x(i)>=0.5 && x(i)<=1.5)
        rho(i) = 1-0.366*(x(i)-0.5);
        t(i) = 1-0.167*(x(i)-0.5);
    elseif (x(i)>=1.5 && x(i)<=3.5)
        rho(i) = 0.634-0.3879*(x(i)-1.5);
        t(i) = 0.833-0.3507*(x(i)-1.5);
    end
end

% Velocity condition

v = 0.59./(rho.*a);

% Initial conditions for Solution vectors
e = t;
p = rho.*t;
u1 = rho.*a;
u2 = rho.*a.*v;
u3 = rho.*a.*((e/(gamma-1))+((gamma/2)*v.^2));

% Time steps
nt = 1400;
%C = input('Input the required Courant no: ');
C = 0.5;
dt = min((C.*dx)./(sqrt(t)+v));

%% Outer time loop
tic
for k = 1:nt
   
    u1_old = u1;
    u2_old = u2;
    u3_old = u3;
    
    % Flux terms in pure form
    
    f1 = u2;
    term1 = (u3./u1) - ((gamma/2)*((u2./u1).^2));
    f2 = ((u2.^2)./u1)+(u1.*term1.*(gamma-1))./(gamma);
    term2 = (gamma*(gamma-1)*(u2.^3))./(2*u1.^2);
    f3 = ((gamma*u2.*u3)./u1)-term2;


    %% Predictor Method
    for i = 2:n-1
        
        % Governing equations
        [du1_dt_p(i),du2_dt_p(i),du3_dt_p(i)] = conservation_p(t,a,i,dx,gamma,f1,f2,f3,rho);
        
        %Solution update
        u1(i) = u1(i) + du1_dt_p(i)*dt;
        u2(i) = u2(i) + du2_dt_p(i)*dt;
        u3(i) = u3(i) + du3_dt_p(i)*dt;
    end
    
    rho = u1./a;
    v = u2./u1;
    t = (gamma-1)*((u3./u1)-((gamma*v.^2)/2));
    
    % Update flux terms in pure form
    
    f1 = u2;
    term1 = (u3./u1) - ((gamma/2)*((u2./u1).^2));
    f2 = ((u2.^2)./u1)+(u1.*term1.*(gamma-1))./(gamma);
    term2 = (gamma*(gamma-1)*(u2.^3))./(2*u1.^2);
    f3 = ((gamma*u2.*u3)./u1)-term2;
    
    % Corrector Method
    for i = 2:n-1
        
        % Governing equations
        [du1_dt_c(i),du2_dt_c(i),du3_dt_c(i)] = conservation_c(t,a,i,dx,gamma,f1,f2,f3,rho);
        
    end
    
    % Compute the average time derivative
    
    du1dt = 0.5*(du1_dt_p + du1_dt_c);
    du2dt = 0.5*(du2_dt_p + du2_dt_c);
    du3dt = 0.5*(du3_dt_p + du3_dt_c);
    
    % Solution update
    for j = 2:n-1
        u1(j) = u1_old(j) + du1dt(j)*dt;
        u2(j) = u2_old(j) + du2dt(j)*dt;
        u3(j) = u3_old(j) + du3dt(j)*dt;
    end
    
    % Apply the boundary conditions
    
    % Inlet B.C.
    u1(1) = rho(1)*a(1);
    u2(1) = 2*u2(2) - u2(3);
    v(1) = u2(1)/u1(1);
    u3(1) = u1(1)*((t(1)/(gamma-1))+((gamma/2)*v(1)^2));
    
    %Outlet B.C.
    u1(n) = 2*u1(n-1) - u1(n-2);
    u2(n) = 2*u2(n-1) - u2(n-2);
    u3(n) = 2*u3(n-1) - u3(n-2);
    
    % Obtaining primitives value i.e. rho, v, t and p
    
    rho = u1./a;
    v = u2./u1;
    t = (gamma-1)*((u3./u1)-((gamma*v.^2)/2));
    p = rho.*t;
    mach = v./(sqrt(t));
    mass = rho.*v.*a;
    
    % Solution at throat (area = 1)
    rho_th(k) = rho(16);
    temp_th(k) = t(16);
    Pr_th(k) = p(16);
    mach_th(k) = mach(16);
    v_th(k) = v(16);
    mass_th(k) = mass(16);
    k_th = 1:nt;
end


% plot(k_th,mass_th,'-b');

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
plot(x,mass); grid on;
xlabel('x/L');
ylabel('{\rho}av/{\rho_o}a_ov_o');
title('Mass flow rate profile');

suptitle('Parameter profile for 1D nozzle flow for Non-Conservation form');

%% Writing output in word file
fid = fopen('conservation.doc ','w');
fprintf(fid,'Parameter for 1D Conservation form nozzle flow \n');
fprintf(fid,'for grid size %d and Courant no %f \n\n', n,C);
fprintf(fid,'Velocity  Density  Temp   Mach no  Massflowrate \n\n');
fprintf(fid, '%f %f %f %f %f \n', [v;rho;t;mach;mass]);
fclose(fid);

%plot(x,mass,'-r');
toc