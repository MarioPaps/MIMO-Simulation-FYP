%Marios Papadopoulos, EE4, 2022, Imperial College.
% 17/01/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the number of co-channel users
% Inputs
% Rxx (FxF Complex) = covariance matrix of the vector of  channel symbol chips received
% N (Integer) = number of antennas used
% L (integer) = number of channel symbols received 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% no_users(scalar)= number of users found 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[no_users]= EstimateNumUsers(Rxx,N,L)
    d= ( eig(Rxx));
    d= abs(sort(d,'descend')); %abs is needed to remove imaginary part from eigenvalue matrix
    LF=zeros(1,N);
    AIC=zeros(1,N);
    MDL=zeros(1,N);
    for k=0:N-1
        prod=1;
        eigen_sum= 0;
        exponent=(1/(N-k));
        for l=k+1:N
            temp=d(l)^(exponent);
            prod=prod*temp;
            eigen_sum=eigen_sum+d(l);
        end
        eigen_sum=(1/(N-k))*eigen_sum;
    
        LF(k+1)=(N-k)*L*(log(prod./eigen_sum)); 
        AIC(k+1)=-2*LF(k+1) + 2*k*(2*N-k);
        MDL(k+1)=-LF(k+1) + 0.5 * k * (2*N-k)*log(L);
    end
    [~,pos_min_AIC]=min(AIC);
    [~,pos_min_MDL]=min(MDL);
    if(pos_min_AIC==pos_min_MDL)
       no_users=pos_min_AIC-1;
    else
        no_users=pos_min_MDL-1;
    end

end