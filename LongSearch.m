function [alphastar,falphastar,gfalphastar,flag]=LS(f,pk,gf,f0,gf0,deltaf,fbar,c1,c2,epsilon)
tau1=9;
alpha0=0;
fi_1=f0;
fprime0=sum(pk.*gf0);
gfi_1=gf0;
alphamax=(fbar-f0)/(c1*fprime0);
%falphamax=f(alphamax);
%fprimealphamax=fprime(alphamax);
alphai_1=alpha0;
alphai=min(1,-1.01*2*deltaf/fprime0);
%fi=f(alphai);
while(1)
    fi=f(alphai);
    if fi<=fbar
        gfi=gf(alphai);
        alphastar=alphai;
        falphastar=fi;
        gfalphastar=gfi;
        flag=0;
        break;
    end
    if fi>f0+c1*alphai*fprime0 || fi>=fi_1
        [alphastar,falphastar,gfalphastar,flag]=zoom(alphai_1, alphai, f, pk, gf, fi_1, fi, gfi_1, f0, fprime0, c1, c2, epsilon);
        break;
    end
    gfi=gf(alphai);
    fprimei=sum(pk.*gfi);
    if abs(fprimei)<=-c2*fprime0
        alphastar=alphai;
        falphastar=fi;
        gfalphastar=gfi;
        flag=0;
        break;
    end
    if fprimei>=0
        [alphastar,falphastar,gfalphastar,flag]=zoom(alphai, alphai_1, f, pk, gf, fi, fi_1, gfi, f0, fprime0, c1, c2, epsilon);
        break;
    end
    if alphamax<=(2*alphai-alphai_1)
        alphai_1=alphai;
        fi_1=fi;
        gfi_1=gfi;
        alphai=alphamax;
    else
        ss=(1/2)*(2*alphai-alphai_1+min(alphamax,alphai+tau1*(alphai-alphai_1)));
        alphai_1=alphai;
        fi_1=fi;
        gfi_1=gfi;
        alphai=ss;
    end
  
end
end

