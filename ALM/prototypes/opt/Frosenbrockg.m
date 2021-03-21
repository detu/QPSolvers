function [f] = Frosenbrockg(x,y);
%
% function [f] = Frosenbrockg(x,y);
%
% This function returns the function value
% of the (two-dimensional) rosenbrock function, given by:
%
%       f(x,y) = 100*(y - x^2)^2 + (1-x)^2 
%
% If x,y are matrices then f is applied componentwise 
% Frosenbrockg is a version of Frosenbrock for n=2 that is convenient
% for contour plots.

f=100*(y-x.^2).^2+(1-x).^2;
