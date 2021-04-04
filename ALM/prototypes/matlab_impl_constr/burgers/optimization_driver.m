%
%  compute the solution of the optimal control problem
%
%   min  0.5 * int_0^T int_0^1 (y(x,t)^2 - y_des(x,t))^2 dx dt
%          + (0.5*omega) * int_0^T int_0^1 u(x,t)^2  dx dt
%
%
%  where y solves the Burgers equation
%
%    y_t(x,t) - nu*y_xx(x,t) + y_x(x,t)*y(x,t) = u(x,t)   x in (0,1), t in (0,1)
%
%                              y(0,t) = y(1,t) = 0        t in (0,1)
%
%                                       y(x,0) = y0(x)    x in (0,1)
%                                      
%
%
% Version June 6, 2008
% Matthias Heinkenschloss
%
%
%      

clear all

set(0, 'defaultaxesfontsize',18,'defaultaxeslinewidth',1.2,...
       'defaultlinelinewidth',0.8,'defaultpatchlinewidth',0.8,...
       'defaulttextfontsize',18);
   
addpath ../optimization


% We use global variables to pass information to the
% functions that compute the solution of Burgers equation, etc.

global BURGERS_GLB



% set problem data
nt        = 80;             % number of time intervals
nx        = 80;             % number of spatial intervals
viscosity = 0.01;           % viscosity
omega     = 0.05;           % control penalty parameter
prob_gen(nx, nt, viscosity, omega);
usr_par   = [];



% perform derivative checks
u = zeros((nx+1)*(nt+1),1);
deriv_check( u, 1, usr_par);

   


% Solve the optimal control problem using Newton's method
fprintf(1,'\n\n Solve the optimal control problem using')

% set initial iterate
u = zeros((nx+1)*(nt+1),1);

optimizer = 1;

if( optimizer == 1)
    fprintf(1,' Newton''s method\n')
    % set options for Newton CG method
    options.iprint = 1; % print level, default = 0
    options.maxit  = 20; %  maximum number of iterations, default = 100
    options.gtol   = 1.e-8; % gradient stopping tolerance, default = 1.e-8
    options.stol   = 1.e-8; % stepsize stopping tolerance, default = 1.e-8

    [ u, iter, iflag ]  = newton_cg( u, options, usr_par);
    fprintf(1,' Newton''s method returned after %d iterations with flag %d \n', iter,iflag )

elseif( optimizer == 2)
    fprintf(1,' limited memory BFGS method\n')
    % set options for limited memory BFGS method
    options.iprint = 1; % print level, default = 0
    options.gtol   = 1.e-8; % gradient stopping tolerance, default = 1.e-8
    options.stol   = 1.e-8; % stepsize stopping tolerance, default = 1.e-8
    options.maxit  = 50; % Maximum number of iterations, default 100
    options.L      = 20;  % LBFGS storage limit, default 20
    [ u, iter, iflag ]  = lbfgs( u, options, usr_par);
    fprintf(1,' LBFGS returned after %d iterations with flag %d \n', iter,iflag )

else
    fprintf(1,' optimizer is %d; must be 1, 2 \n', optimizer)
end

% plot results
clf;
usr_par = xnew( u, 0, usr_par );

figure(1)
% plot the control
mesh((0:BURGERS_GLB.Deltax:1), (0:BURGERS_GLB.Deltat:1), reshape(u,nx+1,nt+1)')
xlabel('x')
ylabel('t')
title('Computed optimal control u_*')
%print -depsc burger_optu

figure(2)
% for plotting purposes augment y by the zero boundary condition
yplot = [ zeros(1,nt+1);
          BURGERS_GLB.y;
          zeros(1,nt+1)];
mesh((0:BURGERS_GLB.Deltax:1), (0:BURGERS_GLB.Deltat:1), yplot')
xlabel('x')
ylabel('t')
title('Corresponding state y(u_*)')
%print -depsc burger_opty

figure(3)
% for plotting purposes augment lambda by the zero boundary condition
lplot = [ zeros(1,nt+1);
          BURGERS_GLB.lambda;
          zeros(1,nt+1)];
mesh((0:BURGERS_GLB.Deltax:1), (BURGERS_GLB.Deltat:BURGERS_GLB.Deltat:1), lplot(:,2:nt+1)')
xlabel('x')
ylabel('t')
title('Corresponding adjoint \lambda(u_*)')
%print -depsc burger_optlambda





