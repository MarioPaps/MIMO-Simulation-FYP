%Generate MAI of a given power based on the NFR levels provided
function [MAI,MAI_nfr] = MAIfromNFR(A,NFR,M,Nosymbs)
    
    %NFR= mean(MAI powers)/ desired user power
    userpower= (1/length(A))*A*ctranspose(A);
    MAI_powers= NFR.*userpower; %the logic is to test performance for increasing MAI levels
    radius= sqrt(MAI_powers);
    MAI_nfr= cell(1,length(NFR));
    bitsMAI= round(rand(M-1,2*Nosymbs));

    for ind=1:length(radius)
        MAI=zeros(M-1,length(A));
        for user=2:M
            MAI(user-1,:)= QPSKMod(bitsMAI(user-1,:),radius(ind),deg2rad(10));
        end
        MAI_nfr{ind}= MAI;
        
    end

    
end