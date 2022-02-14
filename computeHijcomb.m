function[Hijcomb]= computeHijcomb(Hij,J)
    Nc=31;
    Nsc=5;
    N=9;

    J_raised= J^(Nc*Nsc);
    J_raised_t= (J')^(Nc*Nsc);
    
    Hijcomb=[kron(eye(N),J_raised_t)*Hij, kron(eye(N),J_raised)*Hij];

end