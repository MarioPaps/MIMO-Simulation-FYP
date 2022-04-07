%objective: compute hikj: 2NNcNsc x 1
%compute h for a single delay and velocity on a subcarrier and path
function[hikj]= findhVectorised(DOA,l_ik,uk,J,Fj,Nsc,array,userpn)
    Fc=2e10;
    lightvel=3e8;
    Ts=0.1e-6;
    N=max(height(array),width(array));
    Nc=length(userpn);
    del=0:Nc*Nsc-1;
    full_range=(0:1:2*Nc*Nsc-1);
    J_raised= J^l_ik;
    constprod= kron([userpn;zeros(Nc,1)], ones(Nsc,1));

    Sikj=computeManifoldRx(DOA,0,array,Fc,Fj,lightvel);

    diff= 2*pi*Fj*(del-l_ik)*Ts;
    diff=diff';
    temp=[exp(1i*diff);zeros(Nc*Nsc,1)];
    aj_lik= constprod .* temp;

    Fscalar=-(1/lightvel)*(Fc+Fj)*uk;
    Fikj= exp(1i*full_range'*2*pi*Fscalar*Ts);

    par= (J_raised*aj_lik).* Fikj;
    hikj= zeros(2*N*Nc*Nsc, width(Sikj));
    for it=1:width(Sikj)
        hikj(:,it)=kron(Sikj(:,it),par);
    end
   

end
