function[Pcompkjun,hkallj]= gAmpSearch2(gamma,Rxx,lambda_min,H,k,K,M,Nsc,N,Next,J,Fjvec,r,user1pn)
        %we have isolated a path of the desired user and 
        %we search over all subcarriers
        gamma_possible= gamma(k,1:Nsc);
        gamma_possible= abs(gamma_possible);
        gammakj=zeros(1,Nsc);

        %store all the manifold vectors of all j on this specific path
        findpos= ((1:Nsc)-1)*K+k;
        hkallj= H(1:N*Next,findpos);

        %find hkj_hat
        hkj_hat=zeros(N*Next,Nsc);
        for j=1:Nsc
            hkj_hat(:,j)= DoppSTARmanifold(60,140,20,J,Fjvec(j),Nsc,r,user1pn);
        end
        %precompute the products
        products=zeros(N*Next,N*Next,Nsc);
        for j=1:Nsc
             products(:,:,j)= (hkallj(:,j))*(hkallj(:,j))';
        end

        %perform cost function search
        for j=1:Nsc
            ksi=[];
            for iter=1:length(gamma_possible)
                Rkjgamma=Rxx-lambda_min-(gamma_possible(iter)^2*products(:,:,j));
                eigenvalues= eig(Rkjgamma);
                [pos_eval,neg_eval,~,~]= evals(eigenvalues);
                ksi(iter)=(sum(pos_eval+1)) + 10*log10(sum(abs(neg_eval)));
            end
            gammakj=min(ksi(iter))/100000; %this is scaled badly

        end
        Pcompkjun=0;










end