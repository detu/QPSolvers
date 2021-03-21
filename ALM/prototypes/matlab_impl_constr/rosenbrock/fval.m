function [fx] = fval(x, usr_par)
%
% 2D Rosenbrock function
fx = (1-x(1))^2 + 100* ( x(2) - x(1)^2)^2;