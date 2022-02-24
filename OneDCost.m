function[cost]= OneDCost(del_est,vel_est,Fjvec,array,Pn,J,userpn,Nsc,K)
     
      ratio=zeros(Nsc,1);  
      cost=[];
      for k=1:K
          count=1;
          for theta=1:360
                for j=1:Nsc
                    l_ik=del_est(k); %if delay is not zero indexed get rid of 1
                    uk= vel_est(k);
                    hj=DoppSTARmanifold(theta,l_ik,uk,J,Fjvec(j),Nsc,array,userpn);
                    ratio(j)= (ctranspose(hj)*hj)/ (ctranspose(hj)*Pn*hj);
                end
                cost(k,count)= (1/Nsc)* abs( sum(ratio) );
                count=count+1;
          end
      end
      cost=abs(cost);

      %vectorised version
%       theta=(1:0.1:360);
%       for k=1:K
%           count=1;
%           for j=1:Nsc
%                 l_ik=del_est(k); %if delay is not zero indexed get rid of 1
%                 uk= vel_est(k);
%                 hj= findhVectorised(theta,l_ik,uk,J,Fjvec(j),array,userpn);
%                 holder= hj.*conj(hj);
%                 c=sum(holder,1);
%                 for iter=1:length(c)
%                     ratio(j)= c(iter)/ (ctranspose(hj(:,iter))*Pn*hj(:,iter)) ; 
%                 end
%           end
%       end
      %vectorised version findhVectorised(DOA,l_ik,uk,J,Fj,Nsc,array,userpn)
%       theta=(1:0.1:360);
%       holderonj=[];
%       den=[];
%       for k=1:K
%           count=1;
%           l_ik=del_est(k); %if delay is not zero indexed get rid of 1
%           uk= vel_est(k);
%           for j=1:Nsc
%               hj= findhVectorised(theta,l_ik,uk,J,Fjvec(j),Nsc,array,userpn);
%               holderonj=[holderonj,hj];
%           end
%           dotp= holderonj .*conj(holderonj);
%           numer=sum(dotp,1);
%           for iter=1:length(holderonj)
%                 den(iter)= holderonj(:,iter)'* Pn * holderonj(:,iter);
%           end
%          % den2= holderonj'*Pn*holderonj;
%           allratiosj=[];
%           for j=1:Nsc
%                 allratios= numer((j-1)*length(theta)+1: j*length(theta)) ./ den((j-1)*length(theta)+1: j*length(theta));
%                 allratiosj=[allratiosj; allratios];
%           end
%           cost(k,:)= (1/Nsc)* sum(allratiosj,1);
%       end
      
end