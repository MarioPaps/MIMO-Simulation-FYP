function[x]= findXv2(ausers,f,gamma,H,J,M,Nsc,Nc,N,Next,K,n,upper,Pnoise)
    count=1;
    product=zeros(2*N*Nc*Nsc,M*Nsc);
    for i=1:M
        for j=1:Nsc
            row_start= (i-1)*N*Next+1;
            row_end= i*N*Next;
            col_start= (j-1)*K+1;
            col_end= j*K;
            Hij= H(row_start:row_end,col_start:col_end);
            [Hijprev,Hijnext]= computeHijcombv2(Hij,J,Nc,N,Nsc);
            firstmat=[Hijprev,Hij,Hijnext];

            if(n==1)
            currsymb= ausers(j,(i-1)*upper+n);
            nextsymb= ausers(j,(i-1)*upper+n+1);
            secondmat=[zeros(K,1);
                (gamma(:,(i-1)*Nsc+j).*f((j-1)*K+1: j*K,(i-1)*upper+n))* currsymb;
                (gamma(:,(i-1)*Nsc+j).*f((j-1)*K+1: j*K,(i-1)*upper+n+1))* nextsymb];
            elseif(n==upper)
            currsymb= ausers(j,(i-1)*upper+n);
            prevsymb= ausers(j,(i-1)*upper+n-1);
            secondmat=[(gamma(:,(i-1)*Nsc+j).*f((j-1)*K+1: j*K,(i-1)*upper+n-1))* prevsymb;
                     (gamma(:,(i-1)*Nsc+j).*f((j-1)*K+1: j*K,(i-1)*upper+n))* currsymb;
                     zeros(K,1)];
            else
            currsymb= ausers(j,(i-1)*upper+n);
            prevsymb= ausers(j,(i-1)*upper+n-1);
            nextsymb= ausers(j,(i-1)*upper+n+1);
            
            secondmat=[(gamma(:,(i-1)*Nsc+j).*f((j-1)*K+1: j*K,(i-1)*upper+n-1))* prevsymb;
                       (gamma(:,(i-1)*Nsc+j).*f((j-1)*K+1: j*K,(i-1)*upper+n))* currsymb;
                       (gamma(:,(i-1)*Nsc+j).*f((j-1)*K+1: j*K,(i-1)*upper+n+1))* nextsymb];
            end

             %multiplication
            product(:,count)= firstmat * secondmat;
            count=count+1;
        end
    end

        x=sum(product,2);

        disp(n);
end