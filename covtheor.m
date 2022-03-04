%objective: compute Rxx_theor
function[Rxx_theor]= covtheor(H,G,J,N,Nc,Nsc,M,Pnoise)
    user=1;
    K=3;
    Next=2*Nc*Nsc;
    Rdes=zeros(N*Next,N*Next);
    Risi=zeros(N*Next,N*Next);
    Rmai= zeros(N*Next,N*Next);
    for j=1:Nsc
        row_start=(user-1)*N*Next+1;
        row_end= user*N*Next;
        col_start=(j-1)*K+1;
        col_end= j*K;
        term= H(row_start:row_end,col_start:col_end)*G(:,col_start:col_end)*ctranspose(H(row_start:row_end,col_start:col_end));
       

        Rdes= Rdes+term;

        H1jcomb= computeHijcomb(H(row_start:row_end,col_start:col_end),J,Nc,N,Nsc);
        term2=H1jcomb*(kron(eye(2),G(:,col_start:col_end)))*ctranspose(H1jcomb);
        Risi= Risi+term2;
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
           
           check=firstmat*kron(eye(K),G(:,g_start:g_stop))*ctranspose(firstmat);
           Rmai=Rmai+check;

        end
    end
    Rnn= Pnoise*eye(N*Next);
    Rxx_theor=Rdes+Risi+Rmai+Rnn;
end












