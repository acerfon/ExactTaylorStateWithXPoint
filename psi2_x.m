function z=psi2_x(x,y,k,lambda)
z=besselj(1,x*sqrt(lambda^2-k^2)).*cos(k*y)...
    +sqrt(lambda^2-k^2)/2*x.*(besselj(0,x*sqrt(lambda^2-k^2))...
    -besselj(2,x*sqrt(lambda^2-k^2))).*cos(k*y);