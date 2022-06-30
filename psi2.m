function z=psi2(x,y,k,lambda)
z=x.*besselj(1,x*sqrt(lambda^2-k^2)).*cos(k*y);