function[cost]= OneDCost(del_est,vel_est,Fjvec,array,Pn,J,userpn,Nsc,K)
     
      ratio=zeros(Nsc,1);  
      cost=zeros(1,360);
      for k=1:K
          for theta=1:360
                for j=1:Nsc
                    l_ik=del_est(k); %if delay is not zero indexed get rid of 1
                    uk= vel_est(k);
                    hj=DoppSTARmanifold(theta,l_ik,uk,J,Fjvec(j),Nsc,array,userpn);
                    ratio(j)= (ctranspose(hj)*hj)/ (ctranspose(hj)*Pn*hj);
                end
                cost(theta)= (1/Nsc)* abs( sum(ratio) );
          end
      end
    
end