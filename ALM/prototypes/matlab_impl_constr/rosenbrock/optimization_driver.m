% optimization_driver
%
% Use Newton-CG and LBFGS method to solve   min f(x)
% where f is the 2D Rosenbrok function

addpath ../optimization

x = [-0.5;0.5];

% set options for Newton CG Method
options.iprint = 1; % default = 0
options.fid    = 1; % file identifier, default = 1
options.maxit  = 20; %  maximim number of iterations, default = 100
options.gtol   = 1.e-10; % gradient stopping tolerance, default = 1.e-8
options.stol   = 1.e-10; % stepsize stopping tolerance, default = 1.e-8

[ x, iter, iflag ]  = newton_cg( x, options, []);
fprintf(1,' Newton iteration returned after %d iterations with flag %d \n', iter,iflag )
fprintf(1,' and solution x = [%12.5e, %12.5e] \n\n\n', x)
    

x = [-0.5;0.5];

options.maxit    = 30; % Maximum number of iterations, default 100
options.L        = 5;  % LBFGS with storage limited 5, default 20
[ x, iter, iflag ]  = lbfgs( x, options, []);
fprintf(1,' LBFGS returned after %d iterations with flag %d \n', iter,iflag )
fprintf(1,' and solution x = [%12.5e, %12.5e] \n\n\n', x)
