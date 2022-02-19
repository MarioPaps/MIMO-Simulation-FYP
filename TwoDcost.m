function[cost,del_est,uk_est]= TwoDcost(K,Nc,Nsc,Fjvec,akj,Fkj,Pn,J)
    Fc=20e9;
    lightvel=3e8;
    Ts=0.1e-6;
    
    rad_vel_range=0:140;
    del=0:Nc*Nsc-1;
    J_powers=cell(1,length(del));
    for delk=del
        J_powers{delk+1}= J^(delk);
    end

    cost=[];
    ratio=zeros(1,Nsc);
    for delk=del
            for uk=rad_vel_range
                for j=1:Nsc
%                     akj= aikj(:,delk+1,j);
%                     Fcurr=Fkj(:,uk+1,j);
%                     prod= J_powers{1,delk+1}*akj.*Fcurr;
%                     numer=ctranspose(prod)*prod;
                    position= (j-1)*(uk+1)+uk+1;
%                     prod=J_powers{1,delk+1}*akj(:,delk+1,j).*Fkj(:,uk+1,j);
                    prod=J_powers{1,delk+1}*akj(:,delk+1,j).*Fkj(:,position);
                    ratio(j)=(ctranspose(prod)*prod)/ (ctranspose(prod)*Pn*prod);
                    
                end
                cost(delk+1,uk+1)= (1/Nsc)* abs(sum(ratio));
            end
            disp(delk);
    end

    %identify maxima
    [M1,del_ind]= max(cost);
    [~,uk_est]= maxk(M1,K);
    del_est= del_ind(uk_est);
    [del_est, index]= sort(del_est);
    uk_est=uk_est(index);
    del_est=del_est-1; % convert to zero-based indexing
    uk_est=uk_est-1;
end