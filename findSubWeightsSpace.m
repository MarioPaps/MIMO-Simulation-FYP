%compute subspace weights for space-only system
function [space_weights] = findSubWeightsSpace(Rxx,DOA,gamma,array,Fjvec,Fc)
    lightvel=3e8;
    K=length(DOA);
    Nsc=length(Fjvec);
    N=max(size(array));
    %we use the gamma of the first user on j=1
    S=computeManifoldRx(DOA.',0,array,Fc,Fjvec(1),lightvel);
    space_weights=zeros(N,K);

    for k=1:K
        R_unwanted=Rxx- S*gamma(k,1)*S';
        [Evecs,Evals]= eig(R_unwanted);
        [ds,idx]=sort(diag(Evals),'descend');
        E = Evecs(:,idx);
        En= E(:,1:4);
        Pjcomp_unwanted= En*inv((En)'*En)* (En)';
        space_weights(:,k)= Pjcomp_unwanted*S(:,k)*inv(S(:,k)'* Pjcomp_unwanted*S(:,k))*gamma(k,1);
       

    end

   
    
end