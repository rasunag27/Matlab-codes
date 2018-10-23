% Solving steady state 2D heat analysis using different Iterative methods
% The equation is d2u/dx2 + d2u/dy2 = 0;

clear all
close all
clc

% Initialization
Lx = 1; Ly = 1;
nx = 21;
ny = 21;
w = 1.5;
dx = Lx/(nx-1); dy = Ly/(ny-1);
x = linspace(0,Lx,nx);
y = linspace(0,Ly,ny);
[X, Y] = meshgrid(x,y);
tol = 1e-4;
error = 9e9;

% Temp Initialization
T = ones(nx,ny);

%Boundary conditions
T(:,1) = 400; % Left
T(:,end) = 800; % Right
T(1,2:end-1) = 600; % Top
T(end,2:end-1) = 900; % Bottom

% Copying T
Told = T;

display('Three methods to solve 2D Steady heat Conduction');
display('1. Jacobi Iteration');
display('2. Gauss siedel Iteration');
display('3. Over relaxation method');
method = input('Choose any of above three methods for solving: ');

% Jacobi method
iter = 1;
tic
% if method == 1||2||3
while(error>tol)
    for i = 2:nx-1
        for j = 2:ny-1
            if method == 1 
                T(i,j) = jacobi(Told,i,j);
            elseif method == 2
                T(i,j) = gauss(T,i,j);
            elseif method == 3
                T(i,j) = overrelax(T,Told,i,j,w);
            else
                break;
            end
        end
    end
    error = max(max(abs(Told - T)));
    Told = T;
    iter = iter+1;
%   C = contourf(X,Y,T,'edgecolor','none');
%   set(gca,'YDIR','reverse');
%   colormap(jet);
%   clabel(C,'FontSize',12);
%   pause(0.01);

end
toc

%% Plotting
C = contourf(X,Y,T,'edgecolor','none');
set(gca,'YDIR','reverse');
colormap(jet);
clabel(C,'FontSize',12);

% Iteration and Time display
display('-----------------------------------');
fprintf('No of iteration taken is %d\n',iter);
fprintf('Time taken for solving is %f\n',toc);
fprintf('The error obtained is %f\n',error);















