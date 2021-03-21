function s=solvetraspcg(C,c,Del,tol,a,b,d,prec,precini)
%
%
% S. Ulbrich, January 2004
%
% This code comes with no guarantee or warranty of any kind.
%
% function s=solvetraspcg(C,c,Del,tol,a,b,d)
%
% Steihaug-CG-method for the approximate solution of the TR-problem
%
%         min Q(s):=c'*s+1/2*s'*C*s   s.t. ||s||<=Del, a<= d.*s <= b
%
% Input:  C       symmetric nxn-matrix (sparse or dense)
%         c       n-vector (column vector)
%         tol     stopping tolerance: the algorithm stops
%                 if ||c+C s||<=tol*max(1e-12,||c||)
%                 or if TR is left
%                 or if nonpositiv curvature is detected
%         Del     trust region radius
%         a,b,d   bounds and scaling for constraint a<= d.*s<= b
%   prec, precini routines to apply/initialize the preconditioner
%
% Output: s      result after termination
%

del=0.01;
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
  yk=yn;
  gk=gk-alk*Cdk;
  hk=feval(prec,gk);
  gnthn=gk'*hk;
  bek=gnthn/gkthk;
  dk=hk+bek*dk;
  gkthk=gnthn;
 else
  dk=sign(g0'*dk)*dk;
  cc=yk'*yk-Del*Del;
  bb=-2*yk'*dk;
  aa=dk'*dk;
  tau=(-bb+sqrt(bb^2-4*aa*cc))/(2*aa);
  yn=yk-tau*dk;
  msk=(d.*yn>=b);
  if any(msk)
   tau=(d(msk).*yk(msk)-b(msk))./(d(msk).*dk(msk));
   yn=yk-(1-del)*tau*dk;
   stop=1;
  end
  msk=(d.*yn<=a);
  if any(msk)
   tau=(d(msk).*yk(msk)-a(msk))./(d(msk).*dk(msk));
   yn=yk-(1-del)*tau*dk;
   stop=1;
  end
  s=yn;
%  fprintf(1,'Negative curvature in CG-iteration %d\n',it);
  return
 end
end
if (norm(yk)>=Del)
 yk=(Del/norm(yk))*yk;
end
sn=yk;
sp=sn;
sp=(1-del)*min(max(a./d,sp),b./d);
msk=(d.*sn>=b);
if any(msk)
 tau=min(b(msk)./(d(msk).*sn(msk)));
 sn=(1-del)*tau*sn;
end
msk=(d.*sn<=a);
if any(msk)
 tau=min(a(msk)./(d(msk).*sn(msk)));
 sn=(1-del)*tau*sn;
end
predn=-(c'*sn+0.5*sn'*C*sn);
predp=-(c'*sp+0.5*sp'*C*sp);
if predp>=predn
 sn=sp;
 predn=predp;
end
s=sn;
