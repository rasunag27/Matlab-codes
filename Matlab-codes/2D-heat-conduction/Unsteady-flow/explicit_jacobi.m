function sol = explicit_jacobi(Told,alpha_x,alpha_y,i,j)

sol = Told(i,j)+alpha_x*(Told(i+1,j)-2*Told(i,j)+Told(i-1,j))...
                        +alpha_y*(Told(i,j+1)-2*Told(i,j)+Told(i,j-1));
end