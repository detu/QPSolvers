%
%  function [ Hv ] = Hessvec( v, u, usr_par) 
%
%     Purpose:
%
%     Compute the product of the Hessian objective function f(y(u),u)
%     times a vector v, where f(y(u),u) is the objective function 
%     in the optimal control of the unsteady Burger's equation. 
%
%
%     Parameters
%
%     On entry:
%
%     v      Direction v.
%            v((i-1)*(nx+1)+1), ..., v(i*(nx+1))
%            direction at time (i-1)*Deltat, i = 1, ..., nt+1
%
%     u      Control  u.
%            u((i-1)*(nx+1)+1), ..., u(i*(nx+1))
%            controls at time (i-1)*Deltat, i = 1, ..., nt+1
%
%     usr_par user defined parameter. Used to pass problem
%            specific information.
% 
%     On return:
%
%     Hv      Value of Hessian times vector product.
%
%
% Version June 6, 2008
% Matthias Heinkenschloss
%
function [ Hv ] = Hessvec( v, u, usr_par) 

% We use global variables to pass information to the
% functions that compute the solution of Burgers equation, etc.

global BURGERS_GLB

nx     = BURGERS_GLB.nx;
nt     = BURGERS_GLB.nt;
Deltat = BURGERS_GLB.Deltat;
dt2    = Deltat/2;
omega  = BURGERS_GLB.omega;

ny     = nx-1;   % number of y-variables per time step
nu     = nx+1;   % number of u-variables per time step

if( ~BURGERS_GLB.is_state_computed )
    % solve Burgers equation. Store the state as a global variable
    [ BURGERS_GLB.y, iflag ]      = state( u, BURGERS_GLB );
    BURGERS_GLB.is_state_computed = 1;
end
if( ~BURGERS_GLB.is_adjoint_computed )
     % solve the adjoint equation. Store the adjoint as a global variable
     [ BURGERS_GLB.lambda, iflag ]   = adjoint( BURGERS_GLB.y, BURGERS_GLB );
     BURGERS_GLB.is_adjoint_computed = 1;
end

% Solve for w

w = zeros(ny, nt+1);
for i=1:nt
   Mat      = BURGERS_GLB.M + dt2*(BURGERS_GLB.A + Nyp(BURGERS_GLB.y(:,i+1)));
   rhs      = -(-BURGERS_GLB.M + dt2*(BURGERS_GLB.A + Nyp(BURGERS_GLB.y(:,i))))*w(:,i) ...
              +dt2*BURGERS_GLB.B*(v((i-1)*nu+1:i*nu) + v(i*nu+1:(i+1)*nu));
   w(:,i+1) = Mat\rhs;
end



% Solve for p

p = zeros(ny, nt+1);

Mat       = BURGERS_GLB.M + dt2*(BURGERS_GLB.A + Nyp(BURGERS_GLB.y(:,nt+1)));
rhs       = dt2*(BURGERS_GLB.M*w(:,nt+1) + Nyp(w(:,nt+1))'*BURGERS_GLB.lambda(:,nt+1));
p(:,nt+1) = Mat'\rhs;

for i = nt:-1:2
   Mat       = BURGERS_GLB.M + dt2*(BURGERS_GLB.A + Nyp(BURGERS_GLB.y(:,i)));
   rhs       = -(-BURGERS_GLB.M + dt2*(BURGERS_GLB.A + Nyp(BURGERS_GLB.y(:,i))))'*p(:,i+1)...
               + Deltat*BURGERS_GLB.M*w(:,i) ...
               + dt2*(Nyp(w(:,i))'*(BURGERS_GLB.lambda(:,i)+BURGERS_GLB.lambda(:,i+1)));
   p(:,i)    = Mat'\rhs;
end




% Assemble Hessian-vector product

Hv = zeros(size(v));

Hv(1:nu) = (dt2*omega)*(BURGERS_GLB.Q*v(1:nu)) ...
               + dt2*(BURGERS_GLB.B'*p(:,2));

for i = 2:nt
        Hv((i-1)*nu+1:i*nu) = (Deltat*omega)*(BURGERS_GLB.Q*v((i-1)*nu+1:i*nu))...
                              + dt2*(BURGERS_GLB.B'*(p(:,i)+p(:,i+1)));
end
i = nt+1;
Hv((i-1)*nu+1:i*nu) = (dt2*omega)*(BURGERS_GLB.Q*v((i-1)*nu+1:i*nu))...
                          + dt2*(BURGERS_GLB.B'*p(:,i));


% End of hessvec.


