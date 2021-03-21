N=40;
x=[-2*N:2*N]./N;
y=[-N:3*N]./N;
[X,Y]=meshgrid(x,y);
Z=Frosenbrockg(X,Y);
hold off
plotedit on
contour(X,Y,Z,[1:10].^2);
hold on
x0=[-1.9;2];
xopt=trdl(x0,'Frosenbrock',1e-8,1);
hold off
