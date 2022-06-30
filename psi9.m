function z=psi9(x,y,k,lambda)
z=x.*bessely(1,x*sqrt(lambda^2-k^2)).*sin(k*y);