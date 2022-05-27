%this function outputs the SU Doppler-STAR-RAKE and MU Doppler-STAR 
%decorrelating weights
function[w_RAKE,w_dec]= weights4749(H,gamma,K,N,Nc,Nsc)

    %user 1 range is row=1:2*N*Nc*Nsc
%     w_RAKE=zeros(2*N*Nc*Nsc,Nsc);
%     for j=1:Nsc
%        w_RAKE(:,j)= H(1:2*N*Nc*Nsc,(j-1)*K+1:j*K) * gamma(:,j); %RAKE weigthts-equation 47
%     end

    %plan b- compute weights per path,normalise and add
    w_RAKE=zeros(2*N*Nc*Nsc,K);
    for k=1:K
        w_RAKE(:,k)=H(1:2*N*Nc*Nsc,k)*gamma(k);
    end
    w_RAKE= sum(w_RAKE,2,'omitnan');
    
    w_dec=0;

end
