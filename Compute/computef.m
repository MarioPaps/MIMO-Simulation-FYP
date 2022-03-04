%objective is to compute f_ij for all n
function[fi]= computef(in,vels,Fjvec,Fc,Tcs,lightvel,K)
     Nsc=length(Fjvec);
     fi=[];
     for n=1:width(in)
         count=1;
         for j=1:Nsc
             F_vels=-(1/lightvel)*(Fc+Fjvec(j))*vels;
             store= exp(1i*2*n*pi*Fjvec(j)*Tcs) * exp(1i*2*n*pi*(F_vels.')*Tcs);
             start= (j-1)*K+1;
             stop= j*K;
             count=count+1;
             fi(start:stop,n)=store;
         end
     end
     disp('cehc');
end