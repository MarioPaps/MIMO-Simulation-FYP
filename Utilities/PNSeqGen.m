%objective: generate a PN-sequence of length Nc for every user
function[c]= PNSeqGen()
    M=5;
    Nc=31;
    coeffs1=[1 0 0 1 0 1]; %D^5+D^2+1- desired user sequence
    coeffs2=[1 0 1 0 0 1]; %D^5+D^3+1
    coeffs3=[1 0 1 1 0 1]; %D^5+D^3+D^2+D+1
    coeffs4=[1 0 1 1 1 1]; %D^5+D^3+D^2+D+1
    coeffs5=[1 1 0 1 1 1]; %D^5+D^4+D^2+D+1
    
    coeffs=[coeffs1; coeffs2; coeffs3; coeffs4; coeffs5];
    c=zeros(Nc,M);
    
    mseq1= fMSeqGen(coeffs1);
    mseq2= fMSeqGen(coeffs2);
    mseq3= fMSeqGen(coeffs3);
    mseq4= fMSeqGen(coeffs4);
    mseq5= fMSeqGen(coeffs5);
    
    k=(1:Nc);
    gold_seq=zeros(Nc,Nc+2);
    for it=1:numel(k)
        gold_seq(:,it)= fGoldSeq(mseq1,mseq2,k(it));
    end
    gold_seq(:,Nc+1)= mseq1;
    gold_seq(:,Nc+2)= mseq2;
    
    %assign balanced to desired user and offsets to interfering users
    [bal_codes,positions]= fFindBalancedGoldCodes(gold_seq);
    
    c(:,1)= gold_seq(:,positions(1));
    c(:,2)= gold_seq(:,8);

    k=(1:Nc);
    gold_seq=zeros(Nc,Nc+2);
    for it=1:numel(k)
        gold_seq(:,it)= fGoldSeq(mseq4,mseq5,k(it));
    end
    gold_seq(:,Nc+1)= mseq4;
    gold_seq(:,Nc+2)= mseq5;

    [bal_codes2,positions2]= fFindBalancedGoldCodes(gold_seq);

    c(:,3)= gold_seq(:,positions2(1));
    c(:,4)= gold_seq(:,4);
    c(:,5)= gold_seq(:,5);
%     c(:,3)= mseq3;      %gold_seq(:,12);
%     c(:,4)= mseq4; % gold_seq(:,20);
%     c(:,5)= mseq5; %gold_seq(:,27);

    c= -2*c+1;
    
end





