function [xn]=trdlg(x0,fgH,tol,Del0,pausefreq)
%
%
% S. Ulbrich, May 2002
%
% This code comes with no guarantee or warranty of any kind.
%
% function [xn]=trdlg(x0,fgH,tol,pausefreq)
%
% Trust-Region Newton method with Dogleg-type solver for TR-subproblems
%
% Input:  x0      starting point
%         fgH     name of a matlab-function [f,g,H]=fgH(x)
%                 that returns value, gradient and Hessian (dense or sparse)
%                 of the objective function depending on the
%                 number of the given ouput arguments
%         tol     stopping tolerance: the algorithm stops
%                 if ||g(x)||<=tol*min(1,||g(x0)||)
%         Del0    initial trust region radius
%       pausefreq algorithm pauses after pausefreq iterations
%
% Output: xn      result after termination
%

% constants for checking decrease ratio
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
 s=solvetrdl(H,g,Del);
 pred=-(g'*s+0.5*s'*H*s);
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
