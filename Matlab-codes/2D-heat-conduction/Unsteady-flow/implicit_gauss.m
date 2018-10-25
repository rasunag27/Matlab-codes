function sol = implicit_gauss(T,r,i,j,T_prev_dt)
                term1 = (1-2*r)/(1+2*r);
                term2 = r/(1+2*r);
                H = (T(i-1,j)+T(i+1,j));
                V = (T(i,j-1)+T(i,j+1));
                sol = (T_prev_dt(i,j)*term1)+(term2*H)+(term2*V);
end