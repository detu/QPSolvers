function [d,ds]=DDs(x,a,b,g)

ds=abs(g);
d=x-a;
msk=(g<0);
d(msk)=b(msk)-x(msk);
msk=(d>1);
d(msk)=1;
ds(msk)=0;
