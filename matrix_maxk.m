%find the k maxima of a matrix for a surf plot
function[maxima,pos]= matrix_maxk(matrix,num)
    
    pos=zeros(num,2);
    [candidates,~]= maxk(matrix,num);
    candidates= max(candidates);
    candidates_vec= reshape(candidates,[],1);
    maxima= maxk(candidates_vec,num);
    for ind=1:length(maxima)
        [pos(ind,1),pos(ind,2)]= find(matrix==maxima(ind));
    end
    
end