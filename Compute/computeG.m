%the function computes G in equation 25 for all users
function[G]=computeG(gamma,K)
    G=zeros(K,K*width(gamma));
    for iter=1:width(gamma)
        hold= gamma(:,iter).*conj(gamma(:,iter));
        start= (iter-1)*K+1;
        stop=iter*K;
        G(:,start:stop)= diag(hold);
    end

end