clear all
close all
clc

% Constants
L = 0.5;
g = 9.81;
k = 10;
m = 0.1;
c = 0.05;
zeta = c / 2*sqrt(m*k);
omega = sqrt(k/m);

% Initial condition
t = [0 2*pi];
x0 = [pi/2 0];

% Calling ode45

[t, x] = ode45(@(t,x) shm_func(t,x,zeta,omega),t,x0);

%% Plotting
figure(1)
plot(t,x(:,1));
hold on; 
plot(t,x(:,2),'k');
xlabel('Time(s)');
ylabel('Amplitude');
legend('Angle(rad)','Angular speed(rad/s)');

%% Animation
figure(2)
%Origin
O = [0 0];
axis(gca,'equal');
axis([-0.7 0.7 -0.7 0.2]);
grid on;

frame = 1;

% Loop for animation
for i=1:length(t)
    % Mass point
    P = L*[sin(x(i)) -cos(x(i))];

    % Circle in origin
    O_circle = viscircles(O,0.01);
    
    % Pendulum
    pend = line([O(1) P(1)],[O(2) P(2)]);
    % Ball
    ball = viscircles(P,0.05);
    
    % Time for plot update
     pause(0.01);
   
     % Capturing frames
     M(frame) = getframe(gcf);
     frame = frame+1;
    
    % Deleting old pendulum track
    if i<length(t)
        delete(pend);
        delete(ball);
        delete(O_circle);
    end
end

% subplot(6,2,5:12);

% Creating movie file for display
    
movie(M);
v = VideoWriter('Pendulum_Motion_Sunag.avi','Uncompressed AVI');
open(v);
writeVideo(v,M);
close(v);






