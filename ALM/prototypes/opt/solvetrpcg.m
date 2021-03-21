function s=solvetrpcg(C,c,Del,tol,prec,precini)
%
%
% S. Ulbrich, May 2002
%
% This code comes with no guarantee or warranty of any kind.
%
% function s=solvetrpcg(C,c,Del,tol,prec,precini)
%
% Steihaug-PCG-method for the approximate solution of the TR-problem
%
%         min Q(s):=c'*s+1/2*s'*C*s   s.t. ||s||<=Del
%
% Input:  C       symmetric nxn-matrix (sparse or dense)
%         c       n-vector (column vector)
%         tol     stopping tolerance: the algorithm stops
%                 if ||c+C s||<=tol*max(1e-12,||c||)
%                 or if TR is left
%                 or if nonpositiv curvature is detected
%   prec, precini routines to apply/initialize preconditioner
%
% Output: s      result after termination
%

g0=c;
yk=0*c;
gk=g0;
feval(precini,C);
hk=feval(prec,gk);
dk=hk;
gkthk=gk'*hk;
it=0;
sh=[];

% main loop
while norm(gk)>tol*max(1e-12,norm(g0))
 it=it+1;
 Cdk=C*dk;
 dkCdk=dk'*Cdk;
 if (dkCdk>0)
  alk=(gkthk)/(dkCdk);
  yn=yk-alk*dk;
  if isempty(sh)&(norm(yn)>=Del)
   c=yk'*yk-Del*Del;
   b=-2*yk'*dk;
   a=dk'*dk;
   tau=(-b+sqrt(b^2-4*a*c))/(2*a);
   sh=yk-tau*dk;
  end
  yk=yn;
  gk=gk-alk*Cdk;
  hk=feval(prec,gk);
  gnthn=gk'*hk;
  bek=gnthn/gkthk;
  dk=hk+bek*dk;
  gkthk=gnthn;
 else
  dk=sign(g0'*dk)*dk;
  c=yk'*yk-Del*Del;
  b=-2*yk'*dk;
  a=dk'*dk;
  tau=(-b+sqrt(b^2-4*a*c))/(2*a);
  yk=yk-tau*dk;
%  fprintf(1,'Negative curvature in CG-iteration %d\n',it);
  return
 end
end
s=yk;
if (~isempty(sh))
 s=(Del/norm(yk))*yk;
 pred1=-(c'*s+0.5*s'*C*s);
 pred2=-(c'*sh+0.5*sh'*C*sh);
 if (pred2>pred1)
  s=sh;
 end
end
