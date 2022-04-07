function[Hijcomb]= computeHijcomb(Hij,J,Nc,N,Nsc)
   
    J_raised= J^(Nc*Nsc);
    J_raised_t= (J')^(Nc*Nsc);
    
    Hijcomb=[kron(eye(N),J_raised_t)*Hij, kron(eye(N),J_raised)*Hij];

end