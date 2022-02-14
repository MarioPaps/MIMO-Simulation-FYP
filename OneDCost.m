function[var]= OneDCost(del_est,vel_est,Fjvec,array,Pn,J,userpn,Nsc,K)
     
      ratio=zeros(Nsc,1);  
      cost=zeros(1,360);
      for k=1:K
          for theta=1:360
                for j=1:Nsc
                    l_ik=del_est(k);
                    uk= vel_est(k);
                    hj=DoppSTARmanifold(theta,l_ik,uk,J,Fjvec(j),array,userpn);
                    ratio(j)= (ctranspose(hj)*hj)/ (ctranspose(hj)*Pn*hj);
                end
                cost(theta)=abs( sum(ratio)/Nsc );
          end
      end
      var=cost;
end