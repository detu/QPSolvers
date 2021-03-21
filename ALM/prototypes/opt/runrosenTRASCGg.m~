N=40;
x=[-2*N:2*N]./N;
y=[-N:3*N]./N;
[X,Y]=meshgrid(x,y);
Z=Frosenbrockg(X,Y);
hold off
plotedit on
contour(X,Y,Z,[1:10].^2);
hold on
a=[-2;-0.2];
b=[0.6;2.1];
plot([a(1),b(1),b(1),a(1),a(1)],[a(2),a(2),b(2),b(2),a(2)]);
x0=[-1.9;2];
xopt=astrdlg(x0,'Frosenbrock',1e-8,10,a,b,1);
hold off
