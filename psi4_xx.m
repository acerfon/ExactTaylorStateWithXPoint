function z=psi4_xx(x,y,lambda)
z=-lambda./(sqrt(x.^2+y.^2)).*sin(lambda*sqrt(x.^2+y.^2))...
    +lambda*x.^2./(x.^2+y.^2).^(3/2).*sin(lambda*sqrt(x.^2+y.^2))...
    -lambda^2*x.^2./(x.^2+y.^2).*cos(lambda*sqrt(x.^2+y.^2));