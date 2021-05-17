
setenv('BLAS_VERSION','/opt/intel/oneapi/mkl/2021.2.0/lib/intel64/libmkl_core.so');
setenv('LAPACK_VERSION','/opt/intel/oneapi/mkl/2021.2.0/lib/intel64/libmkl_core.so')

H = sparse([5 -2 -1; -2 4 3; -1 3 5]);
f = sparse([2; -35; -47]);

%A   = [];
%b   = [];
%Aeq = [];
%beq = [];
%lb  = [];
%ub  = [];

Aeq = sparse(zeros(3,3));
beq = sparse(zeros(3,1));
lb  = sparse(zeros(3,1));
ub  = sparse([10;10;10]);
x0  = sparse(zeros(3,1));
lambda = sparse(zeros(3,1));

%x = quadprog(H,f,A,b,Aeq,beq,lb,ub) 
tic;
x = qpAL(H,f,Aeq,beq,lb,ub,x0,lambda);
toc;

