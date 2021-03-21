function [g] = grad(x, usr_par)

%
% gradient of 2D Rosenbrock function
% fx = (1-x(1))^2 + 100* ( x(2) - x(1)^2)^2;

g = [ -2*(1-x(1)) - 400*x(1)*( x(2) - x(1)^2);
      200*( x(2) - x(1)^2)];