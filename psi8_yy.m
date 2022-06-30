function z=psi8_yy(x,y,k,lambda)
z=-k^2*x.*besselj(1,x*sqrt(lambda^2-k^2)).*sin(k*y);