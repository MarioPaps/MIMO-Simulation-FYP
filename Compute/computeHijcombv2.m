%the function implements equation 19
function[Hiprev,Hinext]= computeHijcombv2(Hij,J,Nc,N,Nsc)
   
    J_raised= J^(Nc*Nsc); %1st part of equation 19
    J_raised_t= (J')^(Nc*Nsc); %2nd part of equation 19
    Hiprev=kron(speye(N),J_raised_t)*Hij;
    Hinext=kron(speye(N),J_raised)*Hij;

   % if(N>10)
%         Hijcomb=[kron(speye(N),J_raised_t)*Hij, kron(eye(N),J_raised)*Hij]; %equation 19
%     %else
    %    Hijcomb=[kron(eye(N),J_raised_t)*Hij, kron(eye(N),J_raised)*Hij]; %equation 19
   % end
end