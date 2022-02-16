
function[gamma_est]= gPhSearch(x,f,ausers,Pcompkjun,hkallj,psi0,N,Next,Nsc)

    gamma_est=zeros(1,Nsc);
    K=3;
    psi=pi/2;
    guess=zeros(1,Nsc);

    for j=1:Nsc
        currPkjun= Pcompkjun((j-1)*N*Next+1: j*N*Next,(j-1)*N*Next+1: j*N*Next);
        guess(j)= conj(ausers(j,1))* (hkallj(:,j))' *( currPkjun*x(:,1));
        paren=zeros(1,width(x));
        for n=1:width(x)
           paren(n)= f((j-1)*K+1,n)* (hkallj(:,j))' *currPkjun *x(:,n);
        end
        phase= angle(paren-guess(j));
        phase= wrapTo2Pi(phase);
        curly_term= mod((1/width(x))*sum(phase),psi-psi0);
        gamma_est(j)= curly_term+guess(j);
    

    end

    disp('phases found');
end