function [X_min, f_min, iter] = BFGS_LS (f, gf, Hf, x0, Stop_tol, LS_tol, varargin)
iter=1;
deltaf=100;
fbar=0;
c1=10^-4;
c2=0.9;
xk=x0;
f0=f(xk);
gf0=gf(xk);
n=length(x0);
Ck=eye(n);
pk=-Ck*gf0;
if pk==0
    X_min=xk;
    f_min=f0;
    return
end
[alphak,f1,gf1,flag]=LS (@(alpha) f(xk+alpha*pk,varargin{:}), pk, @(alpha) gf(xk+alpha*pk,varargin{:}), f0, gf0, deltaf, fbar, c1, c2, LS_tol);
xkplus1=xk+alphak*pk;
if strcmp(flag,'end of LS')==1
    X_min=xkplus1;
    f_min=f1;
    return
end
iter=iter+1;
deltak=xkplus1-xk;
gammak=gf1-gf0;
Ck=(gammak'*deltak)/(gammak'*gammak)*eye(n);
pk=-Ck*gf0;
%pk=-Ck*gf1;
if pk==0
    X_min=xk;
    f_min=f0;
    return
end
deltaf=f0-f1;%%samira
[alphak,f1,gf1,flag]=LS (@(alpha) f(xk+alpha*pk,varargin{:}), pk, @(alpha) gf(xk+alpha*pk,varargin{:}), f0, gf0, deltaf, fbar, c1, c2, LS_tol);
xkplus1=xk+alphak*pk;
deltaf=f0-f1;
if strcmp(flag,'end of LS')==1
    X_min=xkplus1;
    f_min=f1;
    return
end
%%%%deltaf=f0-f1;%%samira
while (norm(gf1)>Stop_tol)
    iter=iter+1;
    deltak=xkplus1-xk;
    gammak=gf1-gf0;
    a=1/(deltak'*gammak);
    Ck=(eye(n)-a*(deltak*gammak'))*Ck*(eye(n)-a*(gammak*deltak'))+a*(deltak*deltak');
    xk=xkplus1;
    f0=f1;
    gf0=gf1;
    pk=-Ck*gf0;
    if pk==0
        X_min=xk;
        f_min=f0;
        return
    end
    [alphak,f1,gf1,flag]=LS (@(alpha) f(xk+alpha*pk,varargin{:}), pk, @(alpha) gf(xk+alpha*pk,varargin{:}), f0, gf0, deltaf, fbar, c1, c2, LS_tol);
    xkplus1=xk+alphak*pk;
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