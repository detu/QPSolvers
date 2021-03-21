function [f,g,H]=fgHminsurf(z);
global hxq hyq hx hy rect zb msk
%
%
% S. Ulbrich, May 2002
%
% This code comes with no guarantee or warranty of any kind.
%
% function [f,g,H]=fgHminsurf(z)
%
% Based on the triangulation and boundary data defined by the
% global variables the function computes the surface of the
% piecewise linear function with z-values z at the interior
% grid points.
% See runminsurf.m for details.
%
% Input:  z       vector of z-values of the surface at the interior
%                 grid points
%
% Output: f       surface
%         g       gradient
%         H       sparse Hessian
%

zn=zb;
zn(msk)=z;
z=zn;
sqr1=sqrt(hxq*hyq+hxq*(z(rect(:,3))-z(rect(:,1))).^2+...
                       hyq*(z(rect(:,2))-z(rect(:,1))).^2);
sqr2=sqrt(hxq*hyq+hxq*(z(rect(:,2))-z(rect(:,4))).^2+...
                       hyq*(z(rect(:,3))-z(rect(:,4))).^2);
f=0.5*sum(sqr1+sqr2);

N=size(z,1);

if (nargout>1)

g11=zeros(N,1); g12=zeros(N,1); g13=zeros(N,1);
g22=zeros(N,1); g23=zeros(N,1); g24=zeros(N,1);

g11(rect(:,1))=-0.5*(hxq*(z(rect(:,3))-z(rect(:,1)))+...
                hyq*(z(rect(:,2))-z(rect(:,1))))./sqr1;
g24(rect(:,4))=-0.5*(hxq*(z(rect(:,2))-z(rect(:,4)))+...
                hyq*(z(rect(:,3))-z(rect(:,4))))./sqr2;
g13(rect(:,3))=0.5*(hxq*(z(rect(:,3))-z(rect(:,1)))./sqr1);
g23(rect(:,3))=0.5*(hyq*(z(rect(:,3))-z(rect(:,4)))./sqr2);
g12(rect(:,2))=0.5*(hyq*(z(rect(:,2))-z(rect(:,1)))./sqr1);
g22(rect(:,2))=0.5*(hxq*(z(rect(:,2))-z(rect(:,4)))./sqr2);
g=g11+g12+g13+g22+g23+g24;
g=g(msk);
end

if (nargout>2)

r=rect(:,1);
c=rect(:,1);
Gij=(0.5*(hxq+hyq)-2*g11(rect(:,1)).*g11(rect(:,1)))./sqr1;
r=[r;rect(:,1);rect(:,2)];
c=[c;rect(:,2);rect(:,1)];
gij=(-0.5*hyq-2*g11(rect(:,1)).*g12(rect(:,2)))./sqr1;
Gij=[Gij;gij;gij];
r=[r;rect(:,1);rect(:,3)];
c=[c;rect(:,3);rect(:,1)];
gij=(-0.5*hxq-2*g11(rect(:,1)).*g13(rect(:,3)))./sqr1;
Gij=[Gij;gij;gij];
r=[r;rect(:,2)];
c=[c;rect(:,2)];
Gij=[Gij;(0.5*hyq-2*g12(rect(:,2)).*g12(rect(:,2)))./sqr1];
r=[r;rect(:,2);rect(:,3)];
c=[c;rect(:,3);rect(:,2)];
gij=(-2*g12(rect(:,2)).*g13(rect(:,3)))./sqr1;
Gij=[Gij;gij;gij];
r=[r;rect(:,3)];
c=[c;rect(:,3)];
Gij=[Gij;(0.5*hxq-2*g13(rect(:,3)).*g13(rect(:,3)))./sqr1];

r=[r;rect(:,2)];
c=[c;rect(:,2)];
Gij=[Gij;(0.5*hxq-2*g22(rect(:,2)).*g22(rect(:,2)))./sqr2];
r=[r;rect(:,2);rect(:,4)];
c=[c;rect(:,4);rect(:,2)];
gij=(-0.5*hxq-2*g22(rect(:,2)).*g24(rect(:,4)))./sqr2;
Gij=[Gij;gij;gij];
r=[r;rect(:,2);rect(:,3)];
c=[c;rect(:,3);rect(:,2)];
gij=(-2*g22(rect(:,2)).*g23(rect(:,3)))./sqr2;
Gij=[Gij;gij;gij];
r=[r;rect(:,3)];
c=[c;rect(:,3)];
Gij=[Gij;(0.5*hyq-2*g23(rect(:,3)).*g23(rect(:,3)))./sqr2];
r=[r;rect(:,3);rect(:,4)];
c=[c;rect(:,4);rect(:,3)];
gij=(-0.5*hyq-2*g23(rect(:,3)).*g24(rect(:,4)))./sqr2;
Gij=[Gij;gij;gij];
r=[r;rect(:,4)];
c=[c;rect(:,4)];
Gij=[Gij;(0.5*(hxq+hyq)-2*g24(rect(:,4)).*g24(rect(:,4)))./sqr2];
H=sparse(r,c,Gij,N,N);
H=H(msk,msk);
end
