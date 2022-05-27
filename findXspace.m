%compute x[n] for space-only MIMO system
function [x] = findXspace(ausers,gamma,DOA,array,Fjvec,n,upper,M,K,Fc)
    Nsc=min(size(ausers));
    N=max(size(array));
    lightvel=3e8;
    overall_sum=zeros(N,M+Nsc); %accumulation across all paths for every user,subcarrier
    for i=1:M
        for j=1:Nsc
           S_paths=computeManifoldRx(DOA(i,:).',0,array,Fc,Fjvec(j),lightvel);
           product=S_paths.*ausers(j,(i-1)*upper+n).* transpose(gamma(1:K,(j-1)*Nsc+j));
           overall_sum(:,(i-1)*M+j)= sum(product,2);
            
        end
    end
    x=sum(overall_sum,2);
    disp(n);

end