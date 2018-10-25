% Solving steady state 2D heat analysis using Crank-Nicolson Implicit
% method.
% % The equation is d2u/dx2 + d2u/dy2 = 0;

clear all
close all
clc

% Space Initialization 
Lx = 1; Ly = 1;
nx = 11;
ny = 11;
dx = Lx/(nx-1); dy = Ly/(ny-1);

% Time Initialization
time = 0.2;
k = 1.65e-4;             % Thermal Diffusivity
dt = 4;
nt = 100;
r = k*dt/(dx^2);
% k1 = r_x/2;
% k2 = k1;


% Mesh
x = linspace(0,Lx,nx);
y = linspace(0,Ly,ny);
[X, Y] = meshgrid(x,y);

% Limits
tol = 1e-4;
error = 9e9;

% Temp Initialization
T = 400.*ones(nx,ny);

% Boundary conditions
T(:,1) = 400; % Left
T(:,end) = 800; % Right
T(1,2:end-1) = 600; % Top
T(end,2:end-1) = 900; % Bottom

% Copying T
Told = T;


display('Three methods to solve 2D Steady heat Conduction');
display('1. Jacobi Iteration');
display('2. Gauss siedel Iteration');
method = input('Choose any of above three methods for solving: ');

% Implicit method
iter = 1;
T_prev_dt = Told;

for k = 1:nt
    while(error > tol)
        for i = 2:nx-1
            for j = 2:ny-1
                if method == 1
                    T(i,j) = implicit_jacobi(Told,r,i,j,T_prev_dt);
                elseif method == 2
                    T(i,j) = implicit_gauss(T,r,i,j,T_prev_dt);
                else
                    break;
                end
            end
        end  
            T;
            error = max(max(abs(Told - T)));
            Told = T;
            iter = iter+1;
    end
    
    % Plotting
%     C = contourf(X,Y,T,'edgecolor','none');
%     set(gca,'YDIR','reverse');
%     colormap(jet);
%     clabel(C,'FontSize',12);
%     title(sprintf('No of iterations: %d',iter));
%     pause(0.003)
    
    %Looping the previous value
    Told = T;
    T_prev_dt = T;
    error = 9e9;
end

%% Plotting
   C = contourf(X,Y,T,'edgecolor','none');
   set(gca,'YDIR','reverse');
   colormap(jet);
   clabel(C,'FontSize',12);
   title(sprintf('No of iterations: %d',iter));
   pause(0.003)

% % Iteration and Time display
display('2D heat conduction Implicit method');
fprintf('No of iteration taken is %d\n',iter);
fprintf('Time taken for solving is %f\n',toc);











