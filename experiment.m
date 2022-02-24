function[cost]= experiment(del_est,vel_est, theta, Fjvec,array,Pn,J,userpn,Nsc,K)
 
    holderonj=[];
    numer=[];
    den=[];
    cost=zeros(K,length(theta));
    for k=1:K
        l_ik=del_est(k);
        uk=vel_est(k);
        for j=1:Nsc
             hj= findhVectorised(theta,l_ik,uk,J,Fjvec(j),Nsc,array,userpn);
             holderonj=[holderonj;hj];
             %dotp=[dotp; hj.*conj(hj)];
             %numer=[numer; sum(dotp)]
             numer(j,:)=sum(hj.*conj(hj),1);
%              for iter=1:length(hj)
%                    den(j,iter)= hj(:,iter)'*Pn*hj(:,iter); 
%              end
           %  den2= hj'*Pn*hj;
             temp=diag(hj'*Pn*hj);
%              temp=transpose(temp);
             den(j,:)= transpose(temp);
%              result = sum(hj.*(Pn*hj),1);
%              den(j,:)=result;
        end
        for ind=1:length(numer)
            ratio= numer(:,ind)./den(:,ind);
            cost(k,ind)= (1/Nsc)*sum(ratio);
        end




    end


    cost=abs(cost);








end