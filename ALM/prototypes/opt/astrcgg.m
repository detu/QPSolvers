function [xj]=astrcgg(x0,fgH,tol,Del0,a,b,pausefreq)
%
%
% S. Ulbrich, January 2004
%
% This code comes with no guarantee or warranty of any kind.
%
% function [xn]=astrcgg(x0,fgH,tol,Del0,a,b,pausefreq)
%
% Affine-Scaling Trust-Region Newton method with Steihaug-CG solver
% for TR-subproblems
%
% Input:  x0       starting point
%         fgH      name of a matlab-function [f,g,H]=fgH(x)
%                  that returns value, gradient and Hessian (dense or sparse)
%                  of the objective function depending on the
%                  number of the given ouput arguments
%         tol      stopping tolerance: the algorithm stops
%                  if ||g(x)||<=tol*min(1,||g(x0)||)
%         Del0     initial trust region radius
%         a,b      lower and upper bound
%       pausefreq  algorithm pauses after pausefreq steps 
%
% Output: xn      result after termination
%

% constants for check of decrease ratio
del=0.001;
eta1=0.001;
eta2=0.8;

hold on
Del=Del0;
xj=x0;
ab=min(b-a);
% shift xj to the interior if necessary
xj=min(max(xj,a+1e-4*ab),b-1e-4*ab);
[f,g,H]=feval(fgH,xj);
[d,ds]=DDs(xj,a,b,g);
n=size(xj,1);
gh=d.*g;
nmgh0=norm(gh);
nmgh=nmgh0;
it=0;

% main loop
while (norm(gh)>tol*max(1,nmgh0))
 it=it+1;
 dh=sqrt(d);
 Dh=spdiags(dh,0,n,n);
 Ds=spdiags(ds,0,n,n);
 Qt=Dh*H*Dh+Ds;
 gt=dh.*g;
 st=solvetrascg(Qt,gt,Del,min(0.5,norm(gh)),a-xj,b-xj,dh);
 pred=-(gt'*st+0.5*st'*Qt*st)+0.5*st'*(ds.*st);
 s=dh.*st;
 [fn]=feval(fgH,xj+s);
 ared=f-fn;
 if ared>eta1*pred
  xn=xj+s;
  plot([xj(1),xn(1)],[xj(2),xn(2)],'o-')
  if (mod(it,pausefreq)==pausefreq-1)
   fprintf(1,'Press key to continue!\n')
   pause
  end
  xj=xn;
  xj=min(max(xj,a+1e-14*ab),b-1e-14*ab);
  [f,g,H]=feval(fgH,xj);
  [d,ds]=DDs(xj,a,b,g);
  gh=d.*g;
  if ared>eta2*pred
   Del=2*Del;
  end
 else
  Del=min(0.5*Del,norm(st));
 end
 fprintf(1,'it=%3.d   f=%e   ||gh||=%e   Del=%5.3f\n',it,f,norm(gh),Del);
end
fprintf(1,'Successful termination with ||gh||<%e*max(1,||gh0||):\n',tol);
