function [continuity,momentum,energy] = non_conservation_c(v,rho,t,a,j,dx,gamma)

        dvdx = (v(j)-v(j-1))/dx;
        drhodx = (rho(j)-rho(j-1))/dx;
        dtdx = (t(j)-t(j-1))/dx;
        dlog_adx = (log(a(j))-log(a(j-1)))/dx;
        
        %Continuty equation
        continuity = -rho(j)*dvdx-rho(j)*v(j)*dlog_adx - v(j)*drhodx;
        
        %Momentum equation
        momentum = -v(j)*dvdx - (1/gamma)*(dtdx + (t(j)/rho(j))*drhodx);
        
        %Energy equation
        energy = -v(j)*dtdx - (gamma-1)*t(j)*(dvdx + v(j)*dlog_adx);
        
        
        
end
