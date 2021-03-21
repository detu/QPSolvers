function [Hv] = H0vec(v, usr_par)
% 
% compute H0*v, where H0 is the initial BFGS matrix
% H0 is a replacement of the inverse of the Hessian of f!
% Here we use H0 = I
%

Hv = v;
