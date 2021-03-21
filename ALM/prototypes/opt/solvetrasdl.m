function s=solvetrasdl(C,c,Del,a,b,d)
%
%
% S. Ulbrich, January 2004
%
% This code comes with no guarantee or warranty of any kind.
%
% function s=solvetrasdl(C,c,Del,a,b,d)
%
% Dogleg-type-method for the approximate solution of the TR-problem
%
%         min Q(s):=c'*s+1/2*s'*C*s   s.t. ||s||<=Del, a<= d.*s <= b
%
% Input:  C       symmetric nxn-matrix (sparse or dense)
%         c       n-vector (column vector)
%         Del     trust-region radius
%         a,b,d   bounds and scaling for constraint a<= d.*s<= b
%
% Output: s      result after termination
%

del=0.01;
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
msk=(d.*sc>=b);
if any(msk)
 tau=min(b(msk)./(1e-20+d(msk).*sc(msk)));
 sc=(1-del)*tau*sc;
end
msk=(d.*sc<=a);
if any(msk)
 tau=min(a(msk)./(-1e-20+d(msk).*sc(msk)));
 sc=(1-del)*tau*sc;
end
predc=-(c'*sc+0.5*sc'*C*sc);
% Candidate 1: Cauchy point
s=sc;
[R,p]=chol(C);
% Compute the Newton step scaled back to fit in trust region
if p==0
 sn=-((c')/R);
 sn=R\(sn');
 del=min(norm(sn),del);
% Scale Newton step to fit into the trust region
 if norm(sn)>Del
  sn=(Del/norm(sn))*sn;
 end
 sp=sn;
% Candidate 2: Projection of the scaled Newton step onto the box
 sp=(1-del)*min(max(a./d,sp),b./d);

% Candidate 3: Newton step scaled back to fit in box and trust region
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
% Pick best candidate
 predn=-(c'*sn+0.5*sn'*C*sn);
 predp=-(c'*sp+0.5*sp'*C*sp);
 if predp>=predn
  sn=sp;
  predn=predp;
 end
 if predn>=predc
  s=sn;
 end
end
