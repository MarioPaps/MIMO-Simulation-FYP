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
        gamma_range=(0:0.05:1);
        %perform cost function search
        for j=1:Nsc
            ksi=zeros(Nsc,1);
            for iter=1:length(gamma_range)
                Rkjgamma=Rxx-lambda_min-(gamma_range(iter)^2*products(:,:,j));
                eigenvalues= eig(Rkjgamma);
                [pos_eval,neg_eval,~,~]= evals(eigenvalues);
                ksi(iter)=(sum(pos_eval+1)) + 10*log10(sum(abs(neg_eval)));
            end
            gammakj(j)=min(ksi);%/100000; %this is scaled badly
            gammakj(j)= gamma_possible(j);
            
Rkjun((j-1)*N*Next+1: j*N*Next, (j-1)*N*Next+1: j*N*Next)= Rxx-(gammakj(j)^2)*products(:,:,j);
holder=Rkjun((j-1)*N*Next+1: j*N*Next, (j-1)*N*Next+1: j*N*Next);
Pkjun= findPn(holder,M);
Pcompkjun((j-1)*N*Next+1: j*N*Next,(j-1)*N*Next+1: j*N*Next) = eye(N*Next)-Pkjun;
        end
%         Pcompkjun=0;


end