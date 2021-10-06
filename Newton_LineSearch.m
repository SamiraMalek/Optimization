function [X_min, f_min, iter] = Newton_LS (f, gf, Hf, x0, Stop_tol, LS_tol, varargin)
iter=1;
deltaf=100;
fbar=0;
c1=10^-4;
c2=0.9;
xk=x0;
Hk=Hf(xk);
f0=f(xk);
% mineigHk=min(eig(Hk))
% if mineigHk<0
%     Hk=Hk+(-mineigHk+1)*eye(length(Hk));
%     1
% end
gf0=gf(xk);
pk=-Hk^-1*gf0;
if pk==0
    X_min=xk;
    f_min=f0;
    return
end
[alphak,f1,gf1,flag]=LS (@(alpha) f(xk+alpha*pk,varargin{:}), pk, @(alpha) gf(xk+alpha*pk,varargin{:}), f0, gf0, deltaf, fbar, c1, c2, LS_tol);
xkplus1=xk+alphak*pk;
%f1=f(xkplus1);
deltaf=f0-f1;
if strcmp(flag,'end of LS')==1
    X_min=xkplus1;
    f_min=f1;
    return
end
while (norm(xkplus1-xk)>Stop_tol)
    iter=iter+1;
    xk=xkplus1;
    Hk=Hf(xk);
%     mineigHk=min(eig(Hk));
%     if mineigHk<0
%         Hk=Hk+(-mineigHk+1)*eye(length(Hk));
%     end
    f0=f1;
    gf0=gf1;
    pk=-Hk^-1*gf0;
    if pk==0
        X_min=xk;
        f_min=f0;
        return
    end
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
