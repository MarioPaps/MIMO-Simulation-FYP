function [w_j] = findWeights(Rxx,Hj,Gj,gammaj,K,N,Nc,Nsc)
    w_j=zeros(2*N*Nc*Nsc,K);
    for k=1:K
        R_unwanted= Rxx- Hj(:,k)*Gj(k,k)*ctranspose(Hj(:,k));
        if(N>10)
               [Evecs,Evals]= eigs(R_unwanted,50);
        else
               [Evecs,Evals]= eig(R_unwanted);
        end
        [ds,idx]=sort(diag(Evals),'descend');
        E = Evecs(:,idx);
        En= E(:,1:41);

        Pjcomp_unwanted= En*inv((En)'*En)* (En)';
        w_j(:,k)=Pjcomp_unwanted * Hj(:,k)*( (Hj(:,k)'*Pjcomp_unwanted*Hj(:,k)) \ (gammaj(k)));
         %w_j=   Pjcomp_unwanted * Hj*( (Hj'*Pjcomp_unwanted*Hj) \ (gammaj(k)));

    end




end