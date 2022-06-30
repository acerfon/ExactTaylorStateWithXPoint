function z=psi3_y(x,y,k,lambda)
z=-k*x.*bessely(1,x*sqrt(lambda^2-k^2)).*sin(k*y);