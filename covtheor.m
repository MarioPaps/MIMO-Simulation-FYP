%the function computes the theoretical covariance matrix in equation 24
function[Rxx_theor,Rdes,Rmai,Risi,Rnn]= covtheor(H,G,N,Nc,K,Nsc,M,Pnoise)
    user=1;
    Next=2*Nc*Nsc;
    J = [zeros(1,2*Nc*Nsc-1) 0; eye(2*Nc*Nsc-1), zeros(2*Nc*Nsc-1,1)];
    Rdes=zeros(N*Next,N*Next);
    Risi=zeros(N*Next,N*Next);
    Rmai= zeros(N*Next,N*Next);
    for j=1:Nsc
        row_start=(user-1)*N*Next+1;
        row_end= user*N*Next;
        col_start=(j-1)*K+1;
        col_end= j*K;
        term=H(row_start:row_end,col_start:col_end)*G(:,col_start:col_end)*ctranspose(H(row_start:row_end,col_start:col_end)) ;
       
        Rdes= Rdes+term; %desired term

        H1jcomb= computeHijcomb(H(row_start:row_end,col_start:col_end),J,Nc,N,Nsc);
        term2=H1jcomb*(kron(speye(2),G(:,col_start:col_end)))*ctranspose(H1jcomb);
        Risi= Risi+term2; %R_ISI
    end

    for i=2:M
        for j=1:Nsc
           row_start= (i-1)*N*Next+1;
           row_end= i*N*Next;
           col_start= (j-1)*K+1;
           col_end= j*K;
           g_start= col_start+(i-1)*K*Nsc;
           g_stop= col_end+(i-1)*K*Nsc;
            
           Hij=H(row_start:row_end,col_start:col_end);
           firstmat=[computeHijcomb(Hij,J,Nc,N,Nsc),Hij];
           
           check=firstmat*kron(speye(K),G(:,g_start:g_stop))*ctranspose(firstmat);
           Rmai=Rmai+check; %R_MAI

        end
    end
    Rnn= Pnoise*speye(N*Next); %R_nn
    Rxx_theor=Rdes+Risi+Rmai+Rnn; %equation 24
end












