function z=psi8(x,y,k,lambda)
z=x.*besselj(1,x*sqrt(lambda^2-k^2)).*sin(k*y);