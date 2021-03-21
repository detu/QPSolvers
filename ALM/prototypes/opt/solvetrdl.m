function s=solvetrdl(C,c,Del)
%
%
% S. Ulbrich, January 2004
%
% This code comes with no guarantee or warranty of any kind.
%
% function s=solvetrdl(C,c,Del)
%
% Dogleg-type-method for the approximate solution of the TR-problem
%
%         min Q(s):=c'*s+1/2*s'*C*s   s.t. ||s||<=Del
%
% Input:  C       symmetric nxn-matrix (sparse or dense)
%         c       n-vector (column vector)
%         Del     trust-region radius
%
% Output: s      result after termination
%

g=c;
% compute the Cauchy point q(t*sc)=aa*t+0.5*bb*t^2;

nmg=norm(c);
sc=-(1/nmg)*g;
aa=g'*sc;
bb=sc'*(C*sc);
if (bb<=0)
 tau=Del;
else
 tau=min(Del,-aa/bb);
end
sc=tau*sc;
% Candidate 1: Cauchy point
s=sc;
% If C positive definite: Compute the Newton step
[R,p]=chol(C);
% Candidate 1: Compute the Newton step scaled back to fit in trust region
if p==0
 sn=-((c')/R);
 sn=R\(sn');
% Candidate 2: Scaled Newton step to fit into the trust region
% Candidate 3: Dogleg-Point
 if norm(sn)>Del
  sN=(Del/norm(sn))*sn;
  nsn2=sn'*sn;
  nsc2=sc'*sc;
  a=max(1e-20,nsn2+nsc2-2*sn'*sc);
  b=2*sc'*sn-2*nsc2;
  c=min(0,nsc2-Del^2);
  ald=(-b+sqrt(b^2-4*a*c))/(2*a);
  ald=min(1,ald);
  sdl=(1-ald)*sc+ald*sn;
  sn=sN;
  predn=-(c'*sn+0.5*sn'*C*sn);
  preddl=-(c'*sdl+0.5*sdl'*C*sdl);
 else
  predn=-(c'*sn+0.5*sn'*C*sn);
  preddl=predn;
  sdl=sn;  
 end
% Pick best candidate
 if preddl>predn
  s=sdl;
 else
  s=sn;
 end
end
