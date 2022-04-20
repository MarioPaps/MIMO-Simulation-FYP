%the function implements equation 19
function[Hijcomb]= computeHijcomb(Hij,J,Nc,N,Nsc)
   
    J_raised= J^(Nc*Nsc); %1st part of equation 19
    J_raised_t= (J')^(Nc*Nsc); %2nd part of equation 19

   % if(N>10)
        Hijcomb=[kron(eye(N),J_raised_t)*Hij, kron(eye(N),J_raised)*Hij]; %equation 19
    %else
    %    Hijcomb=[kron(eye(N),J_raised_t)*Hij, kron(eye(N),J_raised)*Hij]; %equation 19
   % end
end