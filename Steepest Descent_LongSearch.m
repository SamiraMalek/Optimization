function [X_min, f_min, iter] = SD_LS (f, gf, x0, Stop_tol, LS_tol, varargin)
iter=1;
deltaf=100;
fbar=0;
c1=10^-4;
c2=0.1;
xk=x0;
gf0=gf(xk);
pk=-gf0;
f0=f(xk);
[alphak,f1,gf1,flag]=LS (@(alpha) f(xk+alpha*pk,varargin{:}), pk, @(alpha) gf(xk+alpha*pk,varargin{:}), f0, gf0, deltaf, fbar, c1, c2, LS_tol);
xkplus1=xk+alphak*pk;
%f1=f(xkplus1);
deltaf=f0-f1;
if strcmp(flag,'end of LS')==1
    X_min=xkplus1;
    f_min=f1;
    return
end
xkplus1=xk+alphak*pk;
while (norm(xkplus1-xk)>Stop_tol)
    iter=iter+1;
    xk=xkplus1;
    f0=f1;
    gf0=gf1;
    pk=-gf0;
    [alphak,f1,gf1,flag]=LS (@(alpha) f(xk+alpha*pk,varargin{:}), pk, @(alpha) gf(xk+alpha*pk,varargin{:}), f0, gf0, deltaf, fbar, c1, c2, LS_tol);
    xkplus1=xk+alphak*pk;
    %f1=f(xkplus1);
    deltaf=f0-f1;
    if strcmp(flag,'end of LS')==1
        X_min=xkplus1;
        f_min=f1;
        return
    end
end
X_min=xkplus1;
f_min=f1;
end