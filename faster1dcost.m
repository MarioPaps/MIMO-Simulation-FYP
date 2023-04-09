%%MATLAB function used to evaluate the objective function for DOA estimation
%del_est: vector of estimated delays
%vel_est: vector of estimated radial velocities
%theta: search range
%Fjvec: subcarrier vector
%array: Rx antenna array coordinates
%J: shifting matrix
%userpn: desired user gold code
%Nsc: number of subcarriers
%K: number of paths per user
%cost: matrix obtained from evaluating the cost function of Equation 42
%DOAest: estimated DOAs
function[cost, DOAest]= faster1dcost(del_est,vel_est, theta, Fjvec,array,Pn,J,userpn,Nsc,K)

    %holderonj=[];
    numer=zeros(Nsc,length(theta));
    den=zeros(Nsc,length(theta));
    cost=zeros(K,length(theta));
    for k=1:K
        for j=1:Nsc
             hj= DoppSTARmanifold(theta,del_est(k),vel_est(k),J,Fjvec(j),Nsc,array,userpn);
             %holderonj=[holderonj;hj];
             %dotp=[dotp; hj.*conj(hj)];
             %numer=[numer; sum(dotp)]
             numer(j,:)=sum(hj.*conj(hj),1);
             temp=diag(hj'*Pn*hj);
             den(j,:)= transpose(temp);
        end
       ratio= numer./den;
       cost(k,:)= (1/Nsc)* sum(ratio,1);

    end
    cost=abs(cost);
    DOAest= findMaxofPath(cost);
    DOAest= theta(DOAest);
    
end