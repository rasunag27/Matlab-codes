function [continuity,momentum,energy] = non_conservation_p(v,rho,t,a,j,dx,gamma)

        dvdx = (v(j+1)-v(j))/dx;
        drhodx = (rho(j+1)-rho(j))/dx;
        dtdx = (t(j+1)-t(j))/dx;
        dlog_adx = (log(a(j+1))-log(a(j)))/dx;
        
        %Continuty equation
        continuity = -rho(j)*dvdx-rho(j)*v(j)*dlog_adx - v(j)*drhodx;
        
        %Momentum equation
        momentum = -v(j)*dvdx - (1/gamma)*(dtdx + (t(j)/rho(j))*drhodx);
        
        %Energy equation
        energy = -v(j)*dtdx - (gamma-1)*t(j)*(dvdx + v(j)*dlog_adx);
        
        
        
end
