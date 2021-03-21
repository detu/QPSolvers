function [Hv] = Hessvec(v, x, usr_par)
% Hessian of 2D Rosenbrock function
% fx = (1-x(1))^2 + 100* ( x(2) - x(1)^2)^2;

  
H = [ 2-400*( x(2) - x(1)^2) + 800*x(1)^2,  -400*x(1);
      -400*x(1),   200   ];
  
Hv = H*v;