%the function implements equation 20
function[gammai]= computegamma(beta,DOD,Fjvec,array,K)
    Fc=20e9;
    light_vel=3e8;
    Nsc=length(Fjvec);
    Nbar=16;
    w_bar= ones(Nbar,1);
    holder=zeros(K,1);
    gammai=[];
    for j=1:Nsc
        for path=1:K
            Scurr=computeManifoldTx(DOD(path),0,array,Fc,Fjvec(j),light_vel) ;
            holder(path)= (Scurr)'*w_bar ;
        end
        gammai(:,j)= beta(:,j).* holder;
    end
end