%objective: generate a PN-sequence of length Nc for every user
M=5;
Nc=31;
coeffs1=[1 0 0 1 0 1]; %D^5+D^2+1- desired user sequence
coeffs2=[1 0 1 0 0 1]; %D^5+D^3+1
coeffs3=coeffs2;
coeffs4=coeffs2;
coeffs5=coeffs2;
coeffs=[coeffs1; coeffs2; coeffs3; coeffs4; coeffs5];
c=zeros(Nc,M);

mseq1= fMSeqGen(coeffs1);
mseq2= fMSeqGen(coeffs2);

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
c(:,3)= gold_seq(:,8);
c(:,4)= gold_seq(:,12);
c(:,5)= gold_seq(:,12);


% %generate m sequence
% 
% for it=1:M
%     c(:,it)=fMSeqGen(coeffs(it,:));
% end
%generate gold sequence
disp('check');
c= -2*c+1;
disp('check2');
save('pnseq.mat','c');





