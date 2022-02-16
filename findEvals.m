function[eigenvalues]= findEvals(Mycell)
    eigenvalues=zeros(2790,75);
    Nsc=5;
    count=1;
    for j=1:Nsc
        for iter=1:75
           eigenvalues(:,count)=eig(Mycell{j,iter});
           count=count+1;
    %        [pos_eval,neg_eval,pos_ind,neg_ind]= evals(eigenvalues);
    %        cost(iter)= (1+sum(pos_eval))+ 10*log10(sum(abs(neg_eval)));
    %        disp(iter);
        end
    %     gammakj= min(cost);
    %     Rkjun= Rxx_prac- (gammakj^2)*hkj*ctranspose(hkj);
    end

end