function[S]= computeManifoldRx(theta_ik,phi_ik,r,Fc,Fj,c)
%arg type is used to control the formula depending on units given
        %N= max(size(r));
       %S=zeros(N,length(theta_ik));
        U_ik=[cosd(theta_ik).*cosd(phi_ik),sind(theta_ik).*cosd(phi_ik),zeros(size(theta_ik)).*sind(phi_ik)].';
        Fconst=(1/c)*((Fc+Fj));
        S=exp(-1i*2*pi*Fconst*r*U_ik);

%         for ind=1:length(theta_ik)
%             U_ik=[cosd(theta_ik(ind)).*cosd(phi_ik), sind(theta_ik(ind)).*cosd(phi_ik), sind(phi_ik)]';
%             Fconst=(1/c)*((Fc+Fj));
%             %pos= (theta_ik==theta);
%             S(:,ind)=exp(-1i*2*pi*Fconst*r*U_ik);
%         end
        
end