function sol = overrelax(T,Told,i,j,w)
sol = (1-w)*Told(i,j)+ w*0.25*(T(i+1,j)+T(i-1,j)+ T(i,j+1)+ T(i,j-1));
                          
end
