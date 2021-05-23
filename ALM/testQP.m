setenv('BLAS_VERSION','/opt/intel/oneapi/mkl/2021.2.0/lib/intel64/libmkl_rt.so');
setenv('LAPACK_VERSION','/opt/intel/oneapi/mkl/2021.2.0/lib/intel64/libmkl_rt.so');
TmpLDPath = getenv("LD_LIBRARY_PATH")
H = [5 -2 -1; -2 4 3; -1 3 5];
f = [2; -35; -47];

%A   = [];
%b   = [];
%Aeq = [];
%beq = [];
%lb  = [];
%ub  = [];

Aeq = zeros(3,3);
beq = zeros(3,1);
lb  = zeros(3,1);
ub  = [10;10;10];
x0  = zeros(3,1);
lambda = zeros(3,1);

%x = quadprog(H,f,A,b,Aeq,beq,lb,ub) 
%x = qpAL(H,f,Aeq,beq,lb,ub,x0,lambda);
x = qpAL(H,f);
