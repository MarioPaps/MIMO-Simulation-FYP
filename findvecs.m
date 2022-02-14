function[akj,Fkj]= findvecs(Fjvec,user1pn,Nc,Nsc,Ts)
    rad_vel_range=0:140;
    del=0:Nc*Nsc-1;
    full_range=(0:1:2*Nc*Nsc-1);
    lightvel=3e8;
    Fc=20e9;
    akj=zeros(2*Nc*Nsc,length(del),Nsc);
    Fkj=zeros(2*Nc*Nsc,length(rad_vel_range),Nsc);

    temp1= kron([user1pn;zeros(Nc,1)], ones(Nsc,1));
     for j=1:Nsc
            Fj=Fjvec(j);
             for l_ik=del
                for uk=rad_vel_range
                     diff= 2*pi*Fj*(del-l_ik)*Ts;
                     diff_pad=[diff'; zeros(Nc*Nsc,1)];
                     temp2= exp(1i*diff_pad);
                     prod= temp1.*temp2;
                     akj(:,l_ik+1,j)=prod; %a_j[l_ik] for a specific l_ik
                     
                     Fscalar=-(1/lightvel)*(Fc+Fj)*uk;
                     %I store a 310x1 for every uk on every j
                     Fikj= exp(1i*full_range'*2*pi*Fscalar*Ts);
                     Fkj(:,uk+1,j)=Fikj;
                end
             end



     end


    

end