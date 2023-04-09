%%function used to obtain Doppler-STAR-subspace weights on the %%j-th subcarrier
%Rxx: practical covariance matrix
%Hj: H matrix on the j-th subcarrier
%Gj: G matrix on the j-th subcarrier
%gammaj: gamma vector on the j-th subcarrier
%K: number of paths
%N: number of antennas
%Nc: PN code length
%Nsc: number of subcarriers
%w_j: Doppler-STAR-subspace weight vector for all paths
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
        En= E(:,1:41); %subspace partitioning is done manually because inspection
        %by MLD/AIC
        
        Pjcomp_unwanted= En*inv((En)'*En)* (En)';
        w_j(:,k)=Pjcomp_unwanted * Hj(:,k)*( (Hj(:,k)'*Pjcomp_unwanted*Hj(:,k))\(gammaj(k)));
         %w_j=   Pjcomp_unwanted * Hj*( (Hj'*Pjcomp_unwanted*Hj) \ (gammaj(k)));
    end
end