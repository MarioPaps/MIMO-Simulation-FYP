%Objective: generate an m-sequence of length Nc
function[m_seq]=fMSeqGen(coeffs)
    m= numel(coeffs)-1;
    Nc= 2^(m)-1;
    Qs=zeros(2^(m),m); %matrix to store all Q cycles
    Qs(1,:)=ones(1,m); %initialise first row to 1s
    indices=find(coeffs)-1; %find non-zero coefficients and use zero-based indexing
    pos= nonzeros(indices); %Qi's to xor
    for ind=1:Nc
        Qin=Qs(ind,:); %input at each stage (previous Qs)
        %Quse=[]; for k=1:numel(pos)
        %     Quse=[Quse Qin(pos(k))];
        % end
        Quse=Qin(pos); %pick up the values of Qin at specified indices
        res= mod(nnz(Quse),2); %xor operation equivalent
        Qs(ind+1,:)=[res Qs(ind,1:m-1)];       
    end
    m_seq=Qs(2:end,m); %final output , alternatively 2:end 
    %m_seq= -2*m_seq+1;
end