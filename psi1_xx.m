function z=psi1_xx(x,y,lambda)
z=(1./x-lambda^2*x).*bessely(1,x*lambda)...
    +1/2*lambda./x.*(bessely(0,x*lambda)-bessely(2,x*lambda));