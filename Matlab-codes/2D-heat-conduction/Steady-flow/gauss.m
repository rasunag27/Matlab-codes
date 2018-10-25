function sol = gauss(T,i,j)
 sol = 0.25*(T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1));
                         
end
