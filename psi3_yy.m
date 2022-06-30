function z=psi3_yy(x,y,k,lambda)
z=-k^2*x.*bessely(1,x*sqrt(lambda^2-k^2)).*cos(k*y);