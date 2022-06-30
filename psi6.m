function z=psi6(x,y,lambda)
z=x.*besselj(1,x*lambda).*y;