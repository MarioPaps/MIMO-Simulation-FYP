
function[gammai]= computegamma(beta,DOD,Fjvec,array)
    Fc=20e9;
    light_vel=3e8;
    Nsc=length(Fjvec);
    K=length(DOD);
    Nbar=16;
    w_bar= ones(Nbar,1);
    holder=zeros(K,1);
    gammai=[];
    for j=1:Nsc
        for path=1:K
            Scurr=computeManifoldTx(DOD(path),0,array,Fc,Fjvec(j),light_vel,"hwl") ;
            holder(path)= (Scurr)'*w_bar ;
        end
        gammai(:,j)= beta(:,j).* holder(path);
    end
end