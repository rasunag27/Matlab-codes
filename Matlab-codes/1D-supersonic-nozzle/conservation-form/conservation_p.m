function [continuity,momentum,energy] = conservation_p(t,a,i,dx,gamma,f1,f2,f3,rho)

        dadx = (a(i+1)-a(i))/dx;
        df1dx = (f1(i+1)-f1(i))/dx;
        df2dx = (f2(i+1)-f2(i))/dx;
        df3dx = (f3(i+1)-f3(i))/dx;
        dlog_adx = (log(a(i+1))-log(a(i)))/dx;
        
        j2 = (rho.*t.*dadx)./gamma; 
        
        % Continuity equation
        continuity = -df1dx;
        
        % Momentum equation
        momentum = -df2dx + j2(i);
        
        %Energy equation
        energy = -df3dx;
end
