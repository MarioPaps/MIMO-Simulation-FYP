%compute Rjj
function [Rjj] = spaceRjj(DOAs,array,Fjvec,Fc,M)
    N=max(size(array));
    K=width(DOAs);
    lightvel=3e8;
    intermediate=zeros(N,N,(M-1).^2);
    count=1;
%     Rjj=zeros(N,N);

    for ind=2:M
        for jt=2:M
            %compute S*S' for every path and sum the resulting matrices
            temp=zeros(N,N,K);
            t1=computeManifoldRx(DOAs(ind,:)',0,array,Fc,Fjvec(1),lightvel);
            t2=computeManifoldRx(DOAs(jt,:)',0,array,Fc,Fjvec(1),lightvel);
            for k=1:K
                temp(:,:,k)=t1(:,k)*(t2(:,k))';
            end
            
%           temp=computeManifoldRx(DOAs(ind,:)',0,array,Fc,Fjvec(1),lightvel).*...
%             (computeManifoldRx(DOAs(jt,:)',0,array,Fc,Fjvec(1),lightvel));
           intermediate(:,:,count)=sum(temp,3);
          count=count+1;
        
        end
    end
    Rjj=sum(intermediate,3);


end