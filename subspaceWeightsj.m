%function to compute subspace weights
function [w_j] = subspaceWeightsj(Rxx,Hj,Gj,gammaj,M,Nsc,Nc,N)
    
%     R_unwanted=Rxx- Hj*Gj*ctranspose(Hj);
%     [Evecs,Evals]= eig(R_unwanted);
%     [ds,idx]=sort(diag(Evals),'descend');
%     E = Evecs(:,idx);
%     evals= svd(R_unwanted); %R_unwanted is conjugate symmetric so svd gives basically the same evals in sorted order
    %by inspection, the first 43 eigenvalues are large and the 44th until
    %the end are significantly smaller (of the sorted eigenvalue vector)
    %here, the largest eigenvalues correspond to the noise subspace
    %and the smallest to the signal since we want to remove the effect of
    %the signal
%     En=E(:,1:43);
%     Es=E(:,44:end);
    %Es=E(:,223:end); %for theoretical cov mat
    %     Es=E(:,1:3*M*Nsc-1);
   %En=E(:,1:43);
   %Pjcomp_unwanted= En* inv(En'*En)*En'; these 2 lines produce same output
   %as going through Es
    R_unwanted=Rxx- Hj*Gj*ctranspose(Hj);
    [Evecs,Evals]= eig(R_unwanted);
    [ds,idx]=sort(diag(Evals),'descend');
    E = Evecs(:,idx);
    En= E(:,1:43);

    %Pjcomp_unwanted= eye(2*N*Nc*Nsc)- Es*inv((Es)'*Es)* (Es)';
    %w_j=   Pjcomp_unwanted * Hj* inv(Hj'*Pjcomp_unwanted*Hj)* gammaj;
    %Pjcomp_unwanted= eye(2*N*Nc*Nsc)- Es* (((Es)'*Es) \(Es)') ; 
    Pjcomp_unwanted= En*inv((En)'*En)* (En)';
    w_j=   Pjcomp_unwanted * Hj*( (Hj'*Pjcomp_unwanted*Hj) \ (gammaj));
   % w_j=  Pjcomp_unwanted * Hj /((Hj'*Pjcomp_unwanted*Hj))*gammaj;

    %W sub = P unwanted c*H des est...
% 54 /(H des est'*P unwanted c*H des est)*Gamma des est; % ...
% [2*RaN*Nc x 1]
    

    disp('run');
end