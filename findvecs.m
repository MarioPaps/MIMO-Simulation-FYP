function[akj,Fkj]= findvecs(Fjvec,user1pn,Nc,Nsc,Ts)
    rad_vel_range=1:140;
    del=(0:Nc*Nsc-1);
    full_range=(0:1:2*Nc*Nsc-1);
    lightvel=3e8;
    Fc=20e9;
    akj=zeros(2*Nc*Nsc,length(del),Nsc);
%     Fkj=zeros(2*Nc*Nsc,length(rad_vel_range),Nsc);
    Fkj=zeros(2*Nc*Nsc,Nsc*length(rad_vel_range));
    values=zeros(Nsc*length(del),length(rad_vel_range));

    temp1= kron([user1pn;zeros(Nc,1)], ones(Nsc,1));
    for j=1:Nsc
        Fj=Fjvec(j);
        for l_ik=del
                 diff= 2*pi*Fj*(del-l_ik)*Ts;
                 diff=transpose(diff);
                 temp2=[exp(1i*diff);zeros(Nc*Nsc,1)];
                 prod= temp1.*temp2;
                 akj(:,l_ik+1,j)=prod; %a_j[l_ik] for a specific l_ik      
        end
        for uk=rad_vel_range
                Fscalar= - ((2*pi*(Fc+Fj)*uk*Ts)/lightvel);
                Fstore= exp(1i*full_range'*Fscalar );
              %  position= (j-1)*(uk)+uk;
              Fkj(:,uk,j)=Fstore;
              %  Fkj(:,position)=Fstore;
        end
    end


%      for j=1:Nsc
%             Fj=Fjvec(j);
%              for l_ik=del
%                 for uk=rad_vel_range
%                      diff= 2*pi*Fj*(del-l_ik)*Ts;
%                      diff=diff';
%                      temp2=[exp(1i*diff);zeros(Nc*Nsc,1)];
%                      prod= temp1.*temp2;
%                      akj(:,l_ik+1,j)=prod; %a_j[l_ik] for a specific l_ik
%                      
% %                      Fscalar=(1/lightvel)*(Fc+Fj)*uk*Ts;
%                      %I store a 310x1 for every uk on every j
%                      Fscalar=- ((2*pi*(Fc+Fj)*uk*Ts)/lightvel);
%                      Fstore= exp(1i*full_range'*Fscalar );
%                      position= (j-1)*(uk+1)+uk+1;
% %                      Fkj(:,uk+1,j)=Fstore;
%                      Fkj(:,position)=Fstore;
%                 end
%              end
% 
% 
% 
%      end
end