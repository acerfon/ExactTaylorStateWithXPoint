function z=psi4_y(x,y,lambda)
z=-lambda*y./(sqrt(x.^2+y.^2)).*sin(lambda*sqrt(x.^2+y.^2));