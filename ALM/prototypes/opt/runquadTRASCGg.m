N=50;
x=[-N:N]./N;
y=[-N:N]./N;
[X,Y]=meshgrid(x,y);
Z=Fquadg(X,Y);
hold off
contour(X,Y,Z,[0.05:0.1:1].^2);
hold on
%a=[0.2;0.05];
a=[0.2;-0.2];
b=[1;0.4];
plot([a(1),b(1),b(1),a(1),a(1)],[a(2),a(2),b(2),b(2),a(2)]);
x0=[0.5;0.3];
xopt=astrcgg(x0,'Fquad',1e-8,1,a,b,1);
hold off
