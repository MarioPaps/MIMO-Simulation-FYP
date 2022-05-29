%MUSIC cost function for space-only MIMO
function [musicfn] = MUSIC_cost(Rxx,array,theta,Fjvec,M,Fc)
    [Evecs,Evals]= eig(Rxx);
    [ds,idx]=sort(diag(Evals),'descend');
    E = Evecs(:,idx);
    Es= E(:,1:M);
    En= E(:,M+1:end);
    Ps= Es*inv(ctranspose(Es)*Es)*ctranspose(Es);
    %Pn= En*inv(ctranspose(En)*En)*ctranspose(En);
    Pn= eye(size(Ps))-Ps;

    lightvel=3e8;
    %S= computeManifoldRx(theta.',0,array,Fc,Fjvec(1),lightvel);
    S=spv(array,directions);
    musicfn=diag(S'*Pn*S);
    musicfn= 1./abs(musicfn);


end