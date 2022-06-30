function z=psi8_y(x,y,k,lambda)
z=k*x.*besselj(1,x*sqrt(lambda^2-k^2)).*cos(k*y);