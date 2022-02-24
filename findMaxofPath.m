function[DOA_est]= findMaxofPath(cost1d)
    K=height(cost1d);
    DOA_est=zeros(1,K);
    for k=1:K
        [~,DOA_est(k)]= max(cost1d(k,:));
    end
end