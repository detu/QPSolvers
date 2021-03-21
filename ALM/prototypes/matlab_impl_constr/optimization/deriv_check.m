%
% function [ ]  = deriv_check( x, hess, usr_par)
%
% Check derivatives using finite differences
%
% The user has provide the functions:
%     function [fx] = fval(x, usr_par)         evaluate f(x)
%     function [gx] = grad(x, usr_par)         evaluate grad_f(x)
%     function [Hv] = Hessvec(v, x, usr_par)   evaluate Hess_f(x) * v
%
% and
%     function [usr_par] = xnew( x, iter, usr_par)
%     function [x1x2] = xprod( x1, x2, usr_par).
%     
% xnew is called whenever a new x is generated and it is called
% before any of the three functions fval grad, Hessvec are called 
% with this new x.
%
% xprod evaluates the inner product of x1 and x2.
%
%
% input
%   x       point at which derivatives are checked.
%
%   hess    flag. hess = 0 no Hessian check will be performed.
%
%   usr_par problem specific information. usr_par is not referenced in
%           newton_cg, but passed to the user provided fucntions fval,
%           grad, ....
%
%
% output: none
%                 
%    
%   Matthias Heinkenschloss
%   Department of Computational and Applied Mathematics
%   Rice University
%   June 6, 2008
%
%
function [ ]  = deriv_check( x, hess, usr_par)


fprintf(1,' Derivative checks \n')

fprintf(1,' Gradient check using finite differences (FDs)\n')
fprintf(1,' FD step size    <grad,v>     FD approx.   absolute error \n')

usr_par = xnew( x, 0, usr_par );
f       = fval( x, usr_par );
g       = grad( x, usr_par );
dir     = rand(size(x));
dg      = xprod(g, dir, usr_par);
delta   = 1.e-1;
for d = 1:9
    delta   = delta/10;
    usr_par = xnew( x+delta*dir, 0, usr_par );
    f1      = fval( x+delta*dir, usr_par );

    fprintf(1,' %12.6e  %12.6e  %12.6e  %12.6e  \n', ...
        delta, dg, (f1-f)/delta, abs(dg - (f1-f)/delta) )
end

if( hess )
    % perform Hessian check
    fprintf(1,'\n Hessian check using finite differences (FDs)\n')
    fprintf(1,' FD step size   absolute error \n')

    usr_par = xnew( x, 0, usr_par );
    g       = grad( x, usr_par );
    dir     = rand(size(x));
    delta   = 1.e-1;
    h       = Hessvec( dir, x, usr_par );
    for d = 1:9
        delta   = delta/10;
        usr_par = xnew( x+delta*dir, 0, usr_par );
        g1      = grad( x+delta*dir, usr_par );
        err     = xprod(h - (g1-g)/delta, h - (g1-g)/delta, usr_par);
        fprintf(1,' %12.6e  %12.6e   \n', delta, err )
    end

    % check selfadjointness of Hessian
    fprintf(1,'\n Check if Hessian is selfadjoint \n')

    usr_par = xnew( x, 0, usr_par );
    v1      = rand(size(x));
    v2      = rand(size(x));
    h       = Hessvec( v1, x, usr_par );
    prod1   = xprod(h, v2, usr_par);
    h       = Hessvec( v2, x, usr_par );
    prod2   = xprod(h, v1, usr_par);
    fprintf(1,'<H*v1,v2> = %12.6e, <H*v2,v1> =  %12.6e   \n', prod1, prod2 )
    fprintf(1,'| <H*v1,v2> - <H*v2,v1> | =  %12.6e   \n', abs(prod1 - prod2) )

end
fprintf(1,' End of derivative checks \n\n')
