%Generate MAI of a given power based on the NFR levels provided
function [outputArg1,outputArg2] = MAIfromNFR(A,NFR,M)
    
    %NFR= mean(MAI powers)/ desired user power
    userpower= (1/length(A))*A*ctranspose(A);
    MAI_powers= NFR.*userpower; %the logic is to test performance for increasing MAI levels
    radius= sqrt(MAI_powers);
    MAI=zeros(M-1,length(A));

    for ind=1:length(radius)
        for user=2:M
            bitsMAI= round(rand(1,2*NoSymbs));
            MAI(user-1,:)= QPSKMod(bitsMAI,radius(ind),deg2rad(0));
        end
    end

    %combined output
    ai1=Demux(A,width(A),Nsc);
    

end