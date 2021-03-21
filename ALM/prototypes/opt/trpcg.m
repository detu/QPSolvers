function [xj]=trpcg(x0,fgH,tol,Del0,prec,precini)
%
%
% S. Ulbrich, May 2002
%
% This code comes with no guarantee or warranty of any kind.
%
% function [xn]=trcg(x0,fgH,tol,Del0)
%
% Trust-Region Newton method with Steihaug-CG solver for TR-subproblems
%
% Input:  x0      starting point
%         fgH     name of a matlab-function [f,g,H]=fgH(x)
%                 that returns value, gradient and Hessian (dense or sparse)
%                 of the objective function depending on the
%                 number of the given ouput arguments
%         tol     stopping tolerance: the algorithm stops
%                 if ||g(x)||<=tol*min(1,||g(x0)||)
%         Del0    initial trust region radius
%   prec, precini routines to apply/initialize preconditioner
%
% Output: xn      result after termination
%

% constants for check of decrease ratio
del=0.001;
eta1=0.001;
eta2=0.8;

Del=Del0;
xj=x0;
[f,g,H]=feval(fgH,xj);
nmg0=norm(g);
nmg=nmg0;
it=0;

% main loop
while (norm(g)>tol*max(1,nmg0))
 it=it+1;
 s=solvetrpcg(H,g,Del,min(0.5,norm(g)),prec,precini);
 pred=-(g'*s+0.5*s'*H*s);
 [fn]=feval(fgH,xj+s);
 ared=f-fn;
 if ared>eta1*pred
  xj=xj+s;
  [f,g,H]=feval(fgH,xj);
  if ared>eta2*pred
   Del=2*Del;
  end
 else
  Del=min(0.5*Del,norm(s));
 end
 fprintf(1,'it=%3.d   f=%e   ||g||=%e   Del=%5.3f\n',it,f,norm(g),Del);
end
fprintf(1,'Successful termination with ||g||<%e*max(1,||g0||):\n',tol);
