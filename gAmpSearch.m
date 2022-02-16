%search performed for a single k
%takes 12 minutes to run
function[Pcompkjun,hkallj]= gAmpSearch(gamma,Rxx,H,k,M,Nsc,N,Next)
    gamma_possible= gamma(k,:);
    gamma_possible= abs(gamma_possible);
    evalsRxx=eig(Rxx);
    lambda_min= min(evalsRxx);
    gammakj=zeros(1,Nsc);

    hkallj= [H(1:N*Next,1),H(1:N*Next,4),H(1:N*Next,7),H(1:N*Next,10),H(1:N*Next,13)];
    products=zeros(N*Next,N*Next,Nsc);
    for j=1:Nsc
         products(:,:,j)= (hkallj(:,j))*(hkallj(:,j))';
    end

    Rkjun=zeros(Nsc*N*Next,Nsc*N*Next);
    Pcompkjun=zeros(Nsc*N*Next,Nsc*N*Next);
    for j=1:Nsc
           allpos=[j,j+M,j+2*M,j+3*M,j+4*M];  
           gammarange= gamma_possible(allpos); 
           hkj=hkallj(j);
           ksi=zeros(1,length(gammarange));
           for iter=1:length(gammarange)
               Rkjgamma=Rxx-lambda_min- (gammarange(iter)^2)*products(:,:,j);
               eigenvalues= eig(Rkjgamma);
              [pos_eval,neg_eval,pos_ind,neg_ind]= evals(eigenvalues);
              ksi(iter)= (sum(pos_eval+1))+ 10*log10(sum(abs(neg_eval)));
           end
           gammakj(j)= min(ksi); 
Rkjun((j-1)*N*Next+1: j*N*Next, (j-1)*N*Next+1: j*N*Next)=Rxx-(gammakj(j)^2)*products(:,:,j);
holder=Rkjun((j-1)*N*Next+1: j*N*Next, (j-1)*N*Next+1: j*N*Next);
           Pkjun= findPn(holder,M);
Pcompkjun((j-1)*N*Next+1: j*N*Next,(j-1)*N*Next+1: j*N*Next) = eye(N*Next)-Pkjun;
    end

    disp('amp complete');

end