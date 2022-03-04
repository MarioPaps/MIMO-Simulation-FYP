function[cost,del_est,uk_est]= TwoDcost(K,Nc,Nsc,delays,akj,Fkj,Pn,J) 
    rad_vel_range=1:140;
    del=0:Nc*Nsc-1;
    J_powers=cell(1,length(del));
    for delk=del
        J_powers{delk+1}= J^(delk);
    end
    cost=zeros(length(del),length(rad_vel_range));
    ratio=zeros(1,Nsc);

    

%     for delk=del
%             for uk=rad_vel_range
%                 for j=1:Nsc
% %                     position= (j-1)*(uk)+uk;
%                     J_raised= J_powers{delk+1};
%                     prod=(J_raised*akj(:,delk+1,j)).*Fkj(:,uk,j);
%                   %  prod=(J_raised*akj(:,delk+1,j)).*Fkj(:,position);
%                     ratio(j)=(ctranspose(prod)*prod)/ (ctranspose(prod)*Pn*prod);
%                     
%                 end
%                 cost(delk+1,uk)= (1/Nsc)* abs(sum(ratio));
%             end
%             disp(delk);
%     end
    
    uk_est=zeros(1,K);
    del_est=delays;
    for it=1:length(delays)
        [~,pos]= max(cost(delays(it)+1,:));
        uk_est(it)=pos;
        %cost(delays(it)+1,pos)=  cost(delays(it)+1,pos)+ 10;
    end

% %     %identify maxima
%     [M1,del_ind]= max(cost);
%     [~,uk_est2]= maxk(M1,K);
%     del_est= del_ind(uk_est2);
%     del_est= unique(del_est,'stable');
% 
%     [del_est, index]= sort(del_est);
%     uk_est=uk_est2(index);
    del_est=del_est-1; % convert to zero-based indexing

end