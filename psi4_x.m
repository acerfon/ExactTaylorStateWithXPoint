function z=psi4_x(x,y,lambda)
z=-lambda*x./(sqrt(x.^2+y.^2)).*sin(lambda*sqrt(x.^2+y.^2));