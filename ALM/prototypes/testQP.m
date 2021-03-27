% example of QP problem

x0  = [2; 0];
up  = [inf; inf];
low = [0; 0];
maxit = 100;
[x,histout,costdata,xhist] = gradproj(x0,@objFunction,up,low,maxit);
