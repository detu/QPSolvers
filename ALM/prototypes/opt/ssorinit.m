function [w]=ssorinit(H)

global D DomL

n=size(H,1);
d=spdiags(H,0);
D=sparse(1:n,1:n,d);
om=1.3;
DomL=(D+om*tril(H,-1));
