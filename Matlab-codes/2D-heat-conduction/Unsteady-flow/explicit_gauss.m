function sol = explicit_gauss(T,Told,alpha_x,alpha_y,i,j)

sol = Told(i,j)+alpha_x*(T(i+1,j)-2*T(i,j)+T(i-1,j))...
                        +alpha_y*(T(i,j+1)-2*T(i,j)+T(i,j-1));
end