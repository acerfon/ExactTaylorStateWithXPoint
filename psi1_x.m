function z=psi1_x(x,y,lambda)
z=bessely(1,x*lambda)...
    +1/2*lambda*x.*(bessely(0,x*lambda)-bessely(2,x*lambda));