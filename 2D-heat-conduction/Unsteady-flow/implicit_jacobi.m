function sol = implicit_jacobi(Told,r,i,j,T_prev_dt)
                term1 = (1-2*r)/(1+2*r);
                term2 = r/(1+2*r);
                H = (Told(i-1,j)+Told(i+1,j));
                V = (Told(i,j-1)+Told(i,j+1));
                sol = (T_prev_dt(i,j)*term1)+(term2*H)+(term2*V);
end