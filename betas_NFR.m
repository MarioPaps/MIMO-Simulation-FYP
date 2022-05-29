%output beta for MAI users for given NFR level
function [beta_NFR] = betas_NFR(MAI_powers,M,Nsc,K)
    
    beta_des= 0.1*ones(K,Nsc); %desired user betas
    beta_MAI=zeros(K,(M-1)*Nsc,length(MAI_powers));
    beta_NFR=cell(1,length(MAI_powers));
    beta_mag=linspace(0.2,1,12);

    for ind=1:length(MAI_powers)
        if(ind==1)
            beta_MAI(:,:,ind)= reshape(repelem(beta_des,M-1),K,[]);
            beta_MAI(:,:,ind)=sqrt(MAI_powers(ind));
        else
            beta_MAI(:,:,ind)= beta_mag(ind-1);
            beta_MAI(:,:,ind)=sqrt(MAI_powers(ind));
        end
        beta_NFR{ind}= [beta_des,beta_MAI(:,:,ind)];
    end
    
    %beta_mag= sqrt(MAI_powers);
    %beta_MAI=cell(1,length(MAI_powers));
    %beta_MAI=zeros(K,(M-1)*Nsc,length(MAI_powers));
%     for ind=1:length(MAI_powers)
%         beta_MAI(:,:,ind)=beta_mag(ind).*exp(1i*2*pi*rand(K,(M-1)*Nsc));
%         %beta_MAI{ind}=..
%     end
%     beta_MAI=zeros(K,(M-1)*Nsc);
%     beta_MAI= beta_mag .*exp(1i*2*pi*rand(size(beta_MAI)));

end
% function [beta_NFR] = betas_NFR(MAI_powers,M,Nsc,K)
%     
%     beta_des= 0.01*ones(K,Nsc); %desired user betas
%     beta_MAI=zeros(K,(M-1)*Nsc,length(MAI_powers));
%     beta_NFR=cell(1,length(MAI_powers));
%     beta_mag=linspace(0.2,1,12);
% 
%     for ind=1:length(MAI_powers)
%         if(ind==1)
%             %beta_MAI(:,:,ind)= reshape(repelem(beta_des,M-1),K,[]);
%             beta_MAI(:,:,ind)= kron(beta_des,ones(1,M-1));
%            % beta_MAI(:,:,ind)=sqrt(MAI_powers(ind));
%             %beta_MAI(:,:,ind)= 0.00001;
%         else
%             beta_MAI(:,:,ind)= beta_mag(ind-1);
%             beta_MAI(:,:,ind)=sqrt(MAI_powers(ind));
%         end
%         beta_NFR{ind}= [beta_des,beta_MAI(:,:,ind)];
%     end
%     
%     %beta_mag= sqrt(MAI_powers);
%     %beta_MAI=cell(1,length(MAI_powers));
%     %beta_MAI=zeros(K,(M-1)*Nsc,length(MAI_powers));
% %     for ind=1:length(MAI_powers)
% %         beta_MAI(:,:,ind)=beta_mag(ind).*exp(1i*2*pi*rand(K,(M-1)*Nsc));
% %         %beta_MAI{ind}=..
% %     end
% %     beta_MAI=zeros(K,(M-1)*Nsc);
% %     beta_MAI= beta_mag .*exp(1i*2*pi*rand(size(beta_MAI)));
% 
% end