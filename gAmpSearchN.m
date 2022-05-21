function[Pcompkjun,hkallj]= gAmpSearchN(gamma,Rxx,lambda_min,H,k,K,M,Nsc,N,Next,J,Fjvec,r,user1pn)
        
    gamma_possible= abs(gamma(k,1:Nsc));
    gammapred=zeros(1,Nsc);

    %store all the manifold vectors of all j on this specific path
    findpos= ((1:Nsc)-1)*K+k;
    hkallj= H(1:N*Next,findpos);

    %precompute the products
    products=zeros(N*Next,N*Next,Nsc);
    for j=1:Nsc
         products(:,:,j)= (hkallj(:,j))*ctranspose(hkallj(:,j));
    end

    x=(0:0.1:(max(gamma_possible)));
    x=(0:0.01:0.3);
    Rxx_offset= Rxx-lambda_min;
    %cost fn search
    for j=1:Nsc
        ksi=[];
        for iter=1:length(x)
            Rkjgamma=Rxx_offset- (x(iter).^2)*products(:,:,j);
            eigenvalues= eig(Rkjgamma);
            [pos_eval,neg_eval,~,~]= evals(eigenvalues);
            ksi(iter)=(sum(pos_eval+1)) + 10*log10(sum(abs(neg_eval)));
        end
       ksi=abs(ksi);
       [~,pos]= max(ksi);
       gammapred(j)= x(pos);
    end
    plot(ksi);
    disp('hey');
    
    %use predicted gammas to construct Pcompkjun
    gammapred= gamma_possible;
    for j=1:Nsc
        Rkjun((j-1)*N*Next+1: j*N*Next, (j-1)*N*Next+1: j*N*Next)= Rxx-(gammapred(j)^2)*products(:,:,j);
        holder=Rkjun((j-1)*N*Next+1: j*N*Next, (j-1)*N*Next+1: j*N*Next);
        Pkjun=  holder* inv(holder'*holder)*holder';
        Pcompkjun((j-1)*N*Next+1: j*N*Next,(j-1)*N*Next+1: j*N*Next) = eye(N*Next)-Pkjun;
    end
%     Pcompkjun=0;
end