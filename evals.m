%input: a covariance matrix
function[pos_eval,neg_eval,pos_ind,neg_ind]= evals(eigenvals)
        pos_eval= eigenvals(eigenvals>0);
        neg_eval= eigenvals(eigenvals<0);
        pos_ind= find(eigenvals>0);
        neg_ind= find(eigenvals<0);
end