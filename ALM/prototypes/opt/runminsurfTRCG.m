global hxq hyq hx hy rect zb msk R

% Solves the minimum surface problem
%
%  min int_R sqrt(1+||grad z(x)||^2) dx   s.t. z in H^1(R),
%                                              z=z_b on bnd(R)
% on the rectangle R=[0,Lx]x[0,Ly] for given z_b in H^1(R).
%
% A robust Newton method with sparse direct solver is used.
%
% For the discretization a regular triangulation with piecewise
% linear finite elements is used. The triangulation is obtained
% by halfing the rectangles of a (nx x ny)-grid.

Lx=1;
Ly=1;
nx=80;
ny=80;
hx=Lx/(nx-1);
hy=Ly/(ny-1);
hxq=hx*hx;
hyq=hy*hy;
N=nx*ny;
[X,Y]=meshgrid([0:nx-1]*hx,[0:ny-1]*hy);

% Set the boundary data
Z=(sin(pi*Y)).^2;
zb=reshape(Z',N,1);
z=zb;
% Set mask for interior grid points
v=[1:N];
msk=((mod(v,nx)>1)&(v>nx)&v<=N-nx);
% Generate the triangulation:
% the vertices lie on a regular (nx x ny)-grid
% vertex v(i,j)=(j-1)*nx+i, 1<=i<=nx, 1<=j<=ny
% has coordinates ((i-1)*hx,(j-1)*hy)
% Rectangle k has vertices rect(k,1),...,rect(k,4)  3---4
% and is cut in half to form the two triangles      |\  |
% as shown in the figure                            | \ |
%                                                   |  \|
%                                                   1---2
Nrect=(nx-1)*(ny-1);
rect=zeros(Nrect,4);
rect(:,1)=([1:Nrect]+floor([0:Nrect-1]/(nx-1)))';
rect(:,2)=rect(:,1)+1;
rect(:,3)=rect(:,1)+nx;
rect(:,4)=rect(:,2)+nx;

% Plot initial surface
plotedit off
Z=(reshape(zb,nx,ny))';
surf(X,Y,Z);
pause(0.1);

zopt=trcg(zb(msk),'fgHminsurf',1e-5,10);
zb(msk)=zopt;

% Plot minimal surface
Z=(reshape(zb,nx,ny))';
surf(X,Y,Z);
