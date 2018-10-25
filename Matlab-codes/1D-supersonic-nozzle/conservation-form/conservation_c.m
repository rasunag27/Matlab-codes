function [continuity,momentum,energy] = conservation_c(t,a,i,dx,gamma,f1,f2,f3,rho)

        dadx = (a(i)-a(i-1))/dx;
        df1dx = (f1(i)-f1(i-1))/dx;
        df2dx = (f2(i)-f2(i-1))/dx;
        df3dx = (f3(i)-f3(i-1))/dx;
        dlog_adx = (log(a(i))-log(a(i-1)))/dx;
        
        j2 = (rho.*t.*dadx)./gamma; 
 
        % Continuity equation
        continuity = -df1dx;
        
        % Momentum equation
        momentum = -df2dx + j2(i);
        
        %Energy equation
        energy = -df3dx;
end
