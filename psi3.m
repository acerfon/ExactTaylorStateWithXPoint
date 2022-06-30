function z=psi3(x,y,k,lambda)
z=x.*bessely(1,x*sqrt(lambda^2-k^2)).*cos(k*y);