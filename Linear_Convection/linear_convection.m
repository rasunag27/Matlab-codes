% Basic Linear Convection programme

% Initialization
L = 1;
user = 'Input value of nx: ';  % User input
nx = input(user);
c = 1;
time = 0.4; % Total run time in sec
dt = 0.01;  % Time step
nt = time/dt; 

% Mesh
x = linspace(0,L,nx);
dx = L/(nx-1);

% Array Initialization and Boundary condition
u = ones(1,nx);

% Boundary condition
for i = 1:nx
    if (x(i)>=0.1 && x(i)<=0.3)
        u(i)=2;
    else
        u(i)=1;
    end
end

uold = u;
u_initial = u;

%Time loop - Time marching solution
for j = 1:nt
    % Space loop(BD)    
    for i = 2:nx
        u(i) = uold(i)-(c*dt/dx)*(uold(i)-uold(i-1));
    end
        uold = u;
  
    %Plotting for linear convection
%     plot(j,u_initial,'-ro','linewidth',2,'markeredgecolor','r'); hold on;
%     plot(j,u,'-b'); hold on;
%     axis([0 L 1 2.1])
%     pause(0.05);
   
end

% plot(j,u,'-bo','markeredgecolor','y','markerfacecolor','r');
% xlabel('X');
% ylabel('Velocity');
% legend('Initial Velocity','Time marching velocity','Final Velocity');

% Plot for initial and final figure
figure
plot(x,u_initial,'-ro','linewidth',1.5,'markeredgecolor','r'); hold on;
plot(x,u,'-bo','linewidth',1.5,'markeredgecolor','y','markerfacecolor','r');
xlabel('X');
ylabel('Velocity');
legend('Initial Velocity','Final Velocity');


% Stability is given by Courant number or CFL number
 
C = c*dt/dx % Should be less than 1 to maintain proper stability.






