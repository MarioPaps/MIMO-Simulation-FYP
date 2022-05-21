%estimate amplitude of gamma for one path and one subcarrier of desired
%user
function [gamma_amp_kj] = gampsearch(Rxx,lambda_min,Hj,gamma,Gj,N,Nc,Nsc,Next)
   
   %adapt environment of reference https://spiral.imperial.ac.uk/bitstream/10044/1/946/1/Super-resoultion%20broad%20null%20beamforming.pdf
   % Rmm= Hj*Gj*Hj';
   %     Rmm_cofactor= Rmm;
   %     Rmm_cofactor(2,:)=[];
   %     Rmm_cofactor(:,2)=[];
   %     a0= det(Rmm)/det(Rmm_cofactor)

   %gamma assumed to be 0.1

   %hk_allj= H(1:N*Next,((1:Nsc)-1)*K+k); %store all the manifold vectors of all j on this specific path
   hk_allj=Hj;
   %precompute the products
   products=zeros(N*Next,N*Next,Nsc);
   for j=1:Nsc
         products(:,:,j)= (hk_allj(:,j))*ctranspose(hk_allj(:,j)); %dimensions 2NNcNsc x 2NNcNsc x Nsc
   end

   %cost function search
   x=(0:0.05:1);
   cost=zeros(Nsc,length(x));
   for j=1:Nsc
       for ind=1:length(x)
           Rkj=Rxx-lambda_min- (x(ind)^2*products(:,:,j));
           eigenvalues=eig(Rkj);
           [pos_eval,neg_eval,~,~]= evals(eigenvalues);
           cost(j,ind)= sum(1+pos_eval)+10*log10(sum(abs(neg_eval)));
       end
       cost=abs(cost);
       figure;
       plot(x,cost);
       figure;
       plot(x,1./cost); 
       [~,pos1]= min(cost);
       [~,pos2]= min(1./cost);
   end

   gamma_amp_kj=1;

end