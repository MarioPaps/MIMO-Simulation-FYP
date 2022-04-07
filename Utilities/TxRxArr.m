function[r,r_bar]=TxRxArr(lightvel,Fc)
    lambda=lightvel/Fc;
    rx_bar=[2.56,2.37,1.81,0.98,0,-0.98,-1.81,-2.37,-2.56,-2.37,-1.81,-0.98,0,0.98,1.81,2.37];
    ry_bar=[0,0.98,1.81,2.37,2.56,2.37,1.81,0.98,0,-0.98,-1.81,-2.37,-2.56,-2.37,-1.81,-0.98];
    rz_bar=zeros(1,length(rx_bar));
    r_bar=[rx_bar;ry_bar;rz_bar];
    %if r is in [m], I divide by (lambda/2) to give it (lambda/2) units
    %if r is in [lambda/2], I multiply by (lambda/2) to give it [m] units
    r_bar=r_bar*(lambda/2); %convert to meters
    
    rx=[1.43 0.90 -0.05 -0.98 -1.45 -1.24 -0.45 0.55 1.29];
    ry=[0.30 1.15 1.46 1.09 0.20 -0.77 -1.39 -1.36 -0.69];
    rz=zeros(1,length(rx));
    r=[rx;ry;rz];
    r= r*(lambda/2);
   
end