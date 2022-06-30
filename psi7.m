function z=psi7(x,y,lambda)
z=x.*bessely(1,x*lambda).*y;