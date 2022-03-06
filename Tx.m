function[Tx_out]=Tx(in,inpn,Nsc,N_bar,Tc)
      Nc=length(inpn);
      Fj=((0:Nsc-1)*(1/Tc))';
      midx=(1:Nc*Nsc);
      dim=max(midx);
      F=zeros(Nsc,dim);
      for midx=1:Nc*Nsc
             tempo= exp(1i*2*pi*Fj*midx*(Tc/Nsc));
             F(:,midx)=tempo;
      end
        
     norms= sqrt(sum(abs(in.^2),1)); %find norms of all columns
     in= bsxfun(@rdivide,in,norms); % Normalize the matrix columns

     B=kron(inpn,ones(Nsc,1))';
     M=cell(1,width(in));
     Fsignal=cell(1,width(in));
     W_bar= ones(N_bar,Nsc);
     for n=1:width(in)
         t= (n*Nc*Nsc+(1:Nc*Nsc))* (Tc/Nsc);
         F_n=[];
         for iter=1: Nc*Nsc
            F_n(:,iter)= exp(1i*2*pi*Fj'*t(iter));
         end
          M{n}= kron(in(:,n),B) .* F_n;
          Fsignal{n}=W_bar*M{n};
     end
     Tx_out= cell2mat(Fsignal);
end