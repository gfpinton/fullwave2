function [a b] = ab2(dx,kappax,alphax,dT)
b=exp(-(dx./kappax+alphax).*dT);
a=dx./(kappax*(dx+kappax*alphax)).*(b-1);

