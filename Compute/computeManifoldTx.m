function[S]= computeManifoldTx(theta_ik,phi_ik,r,Fc,Fj,c)
%arg type is used to control the formula depending on units given
        Nbar= max(size(r));
        S=zeros(Nbar,length(theta_ik));
        for theta=theta_ik
            U_ik=[cosd(theta).*cosd(phi_ik), sind(theta).*cosd(phi_ik), sind(phi_ik)]';
            Fconst=(1/c)*((Fc+Fj));
            pos= (theta_ik==theta);
            S(:,pos)=exp(1i*2*pi*Fconst*r*U_ik);
        end
        
end