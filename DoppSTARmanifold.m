%objective: compute hikj: 2NNcNsc x 1
%compute h for a single delay and velocity on a subcarrier and path
function[hikj]= DoppSTARmanifold(DOA,l_ik,uk,J,Fj,array,userpn)
    Fc=2e10;
    lightvel=3e8;
    Ts=0.1e-6;
    Nsc=5;
    Nc=length(userpn);
    del=0:Nc*Nsc-1;
    full_range=(0:1:2*Nc*Nsc-1);
    J_raised= J^l_ik;
    constprod= kron([userpn;zeros(Nc,1)], ones(Nsc,1));

    Sikj=computeManifoldRx(DOA,0,array,Fc,Fj,lightvel,"hwl");

    diff= 2*pi*Fj*(del-l_ik)*Ts;
    diff_pad=[diff'; zeros(Nc*Nsc,1)];
    temp= exp(1i*diff_pad);
    aj_lik= constprod .* temp;

    Fscalar=-(1/lightvel)*(Fc+Fj)*uk;
    Fikj= exp(1i*full_range'*2*pi*Fscalar*Ts);

    par= J_raised*aj_lik .* Fikj;
    hikj= kron(Sikj,par);

end
