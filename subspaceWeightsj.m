%function to compute subspace weights
function [w_j] = subspaceWeightsj(Rxx,Hj,Gj,gammaj,M,Nsc,Nc,N)
    
    R_unwanted=Rxx- Hj*Gj*ctranspose(Hj);
    [Evecs,Evals]= eig(R_unwanted);
    [ds,idx]=sort(diag(Evals),'descend');
    E = Evecs(:,idx);
    %Es= E(:,end-44);
    %Es= E(:,1:1800);
   % Es= E(:, M*Nsc+1:end); %eigenvectors span the unwanted signal subspace
    Es= E(:,44:end);
    %Es= E(:,1:44);

    %Pjcomp_unwanted= eye(2*N*Nc*Nsc)- Es*inv((Es)'*Es)* (Es)';
    %w_j=   Pjcomp_unwanted * Hj* inv(Hj'*Pjcomp_unwanted*Hj)* gammaj;
    Pjcomp_unwanted= eye(2*N*Nc*Nsc)- Es* (((Es)'*Es) \(Es)') ; 
    w_j=   Pjcomp_unwanted * Hj*( (Hj'*Pjcomp_unwanted*Hj) \ (gammaj));
    

    disp('run');
%     w_j=1;

end