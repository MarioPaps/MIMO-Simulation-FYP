%estimate gamma amplitude on particular path- the function assumes the
%system uses Nsc=1 and that gamma coefficients are real
%pass H in as the matrix of the first user on all paths
function[gamma_amp_est]=gammaAmpSearch(Rxx,gamma,H,k,N,Nc,Nsc)
    Rxx_offset=Rxx-min(eig(Rxx));
    products=zeros(N*2*Nc*Nsc,N*2*Nc*Nsc);
    hkallj= H(:,k);
    products= (hkallj(:,1))*ctranspose(hkallj(:,1));
     
%     for j=1:Nsc
%          products(:,:,j)= (hkallj(:,j))*ctranspose(hkallj(:,j));
%     end

    options = optimset('TolX',1e-3);
    lim=fix(abs(gamma(k,Nsc))*100)/100 + 0.01 ;
    gamma_amp_est=fminbnd(@ksisearch,0,lim,options);
    function ksi= ksisearch(x)
            Rkjgamma=Rxx_offset- (x.^2).*products(:,:,1);
            eigenvalues= eig(Rkjgamma);
            [pos_eval,neg_eval,~,~]= evals(eigenvalues);
            ksi=(sum(pos_eval+1)) + 10*log10(sum(abs(neg_eval)));      
    end
end