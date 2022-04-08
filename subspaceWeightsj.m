%function to compute subspace weights
function [w_j] = subspaceWeightsj(Rxx,Hj,Gj,gammaj,M,Nsc,Nc,N)
    
    R_unwanted=Rxx- Hj*Gj*ctranspose(Hj);
    [Evecs,Evals]= eig(R_unwanted);
    [ds,idx]=sort(diag(Evals),'descend');
    E = Evecs(:,idx);
    %by inspection, the first 43 eigenvalues are large and the 44th until
    %the end are significantly smaller (of the sorted eigenvalue vector)
    %here, the largest eigenvalues correspond to the noise subspace
    %and the smallest to the signal since we want to remove the effect of
    %the signal
    Es=E(:,44:end);
%     Es=E(:,1:3*M*Nsc-1);
    

    %Pjcomp_unwanted= eye(2*N*Nc*Nsc)- Es*inv((Es)'*Es)* (Es)';
    %w_j=   Pjcomp_unwanted * Hj* inv(Hj'*Pjcomp_unwanted*Hj)* gammaj;
    Pjcomp_unwanted= eye(2*N*Nc*Nsc)- Es* (((Es)'*Es) \(Es)') ; 
    w_j=   Pjcomp_unwanted * Hj*( (Hj'*Pjcomp_unwanted*Hj) \ (gammaj));
    

    disp('run');
end