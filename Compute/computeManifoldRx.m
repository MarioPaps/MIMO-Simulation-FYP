function[S]= computeManifoldRx(theta_ik,phi_ik,r,Fc,Fj,c, arg_type)
%arg type is used to control the formula depending on units given
        S=zeros(height(r),length(theta_ik));
        if(arg_type=="hwl")
            for ind=1:length(theta_ik)
                U_ik=[cosd(theta_ik(ind)).*cosd(phi_ik), sind(theta_ik(ind)).*cosd(phi_ik), sind(phi_ik)]';
                Fconst=(1/c)*((Fc+Fj));
                %pos= (theta_ik==theta);
                S(:,ind)=exp(-1i*2*pi*Fconst*r*U_ik);
            end
        end
end