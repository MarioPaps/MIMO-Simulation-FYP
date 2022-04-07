function[ds_mag]= manifoldDer(theta_ik,phi_ik,array,Fc,Fj,c, arg_type)
      
     if(arg_type=="hwl")
         for ind=1:length(theta_ik)
           U_ik=[cosd(theta_ik(ind)).*cosd(phi_ik), sind(theta_ik(ind)).*cosd(phi_ik), sind(phi_ik)]';
           derivative=[-sind(theta_ik); cosd(theta_ik);0];
           store=(array*derivative);
           main= exp(-1i*2*pi*((Fc+Fj)/c)*array*U_ik);
           scalar=-1i*2*pi*((Fc+Fj)/c);
           %ds= exp(-1i*2*pi*((Fc+Fj)/c)*array*U_ik)*(-1i*2*pi*((Fc+Fj)/c)).*(array*derivative);
           ds= main.*store;
           ds_mag= norm(ds);
         end
     end

end