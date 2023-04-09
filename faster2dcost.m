%%MATLAB function to implement the 2-D delay-radial velocity cost function
%Nc: PN code length
%Nsc: number of subcarriers
%K: number of paths per user
%akj: vector defined in Equation 26
%Fkj: vector defined in Equation 27
%Pn: projection operator onto noise subspace
%J: shifting matrix
%cost: matrix obtained from evaluating the cost function of Equation 42
function[cost,del_est,uk_est]= faster2dcost(K,Nc,Nsc,delays,rad_vel_range,akj,Fkj,Pn,J) 
    
    del=0:Nc*Nsc-1;
    J_powers=cell(1,length(del));
    for delk=del
        J_powers{delk+1}= J^(delk);
    end
    cost=zeros(length(del),length(rad_vel_range));
    
    for delk=del
        ratio=zeros(Nsc,length(rad_vel_range));
        for j=1:Nsc
            tempo= J_powers{delk+1} * akj(:,delk+1,j);
            tempor= tempo .*Fkj(:,:,j);
            numer= sum(tempor.*conj(tempor),1);
            den= diag(tempor'*Pn*tempor);
            den= transpose(den);
            ratio(j,:)=numer./den;
        end
        cost(delk+1,:)= (1/Nsc)*sum(ratio,1);

    end
    cost=abs(cost);
    uk_est=zeros(1,K);
    del_est=delays;
    for it=1:length(delays)
        [~,pos]= max(cost(delays(it)+1,:));
        uk_est(it)=rad_vel_range(pos);
        %cost(delays(it)+1,pos)=  cost(delays(it)+1,pos)+ 10;
    end
end