function sol = jacobi(Told,i,j)

    sol = 0.25*(Told(i+1,j)+Told(i-1,j)+Told(i,j+1)+Told(i,j-1));
    
end