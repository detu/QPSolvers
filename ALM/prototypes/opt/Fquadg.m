function f=Fquadg(x,y)
%
% function f=Fquadg(x,y)
%
% This function returns the function value
% of the (two-dimensional) quadratic function
%
%       f(x,y) = 0.5*(x^2+20*y^2)
%
% If x,y are matrices then f is applied componentwise.
%
% Fquadg is a version of Fquad for n=2 that is convenient
% for contour plots.

f=0.5*(x.^2+20*y.^2);
