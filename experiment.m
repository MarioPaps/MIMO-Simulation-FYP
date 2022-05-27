function[cost, DOAest]= experiment(del_est,vel_est, theta, Fjvec,array,Pn,J,userpn,Nsc,K)

    holderonj=[];
    numer=zeros(Nsc,length(theta));
    den=zeros(Nsc,length(theta));
    cost=zeros(K,length(theta));
    for k=1:K
        for j=1:Nsc
             hj= DoppSTARmanifold(theta,del_est(k),vel_est(k),J,Fjvec(j),Nsc,array,userpn);
             holderonj=[holderonj;hj];
             %dotp=[dotp; hj.*conj(hj)];
             %numer=[numer; sum(dotp)]
             numer(j,:)=sum(hj.*conj(hj),1);
%              for iter=1:length(hj)
%                    den(j,iter)= hj(:,iter)'*Pn*hj(:,iter); 
%              end
           %  den2= hj'*Pn*hj;
             temp=diag(hj'*Pn*hj);
             den(j,:)= transpose(temp);
        end
       ratio= numer./den;
       cost(k,:)= (1/Nsc)* sum(ratio,1);

    end
    cost=abs(cost);
    DOAest= findMaxofPath(cost);
    DOAest= theta(DOAest);
    
end