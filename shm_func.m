function [xdot] = shm_func(t,x,zeta,omega)

xdot = zeros(2,1);

xdot(1) = x(2);

xdot(2) = -2*zeta*omega*x(2) - omega^2*x(1);

end