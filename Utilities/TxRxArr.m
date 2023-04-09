%UCA:string parameter to generate two different Rx arrays
function[r,r_bar]=TxRxArr(lightvel,Fc,UCA)
    lambda=lightvel/Fc;
    rx_bar=[2.56,2.37,1.81,0.98,0,-0.98,-1.81,-2.37,-2.56,-2.37,-1.81,-0.98,0,0.98,1.81,2.37]';
    ry_bar=[0,0.98,1.81,2.37,2.56,2.37,1.81,0.98,0,-0.98,-1.81,-2.37,-2.56,-2.37,-1.81,-0.98]';
    rz_bar=zeros(length(rx_bar),1);
    r_bar=[rx_bar,ry_bar,rz_bar];
    %if r is in [m], I divide by (lambda/2) to give it (lambda/2) units
    %if r is in [lambda/2], I multiply by (lambda/2) to give it [m] units
    r_bar=r_bar*(lambda/2); %convert to meters

    if(UCA==9)
        rx=[1.43 0.90 -0.05 -0.98 -1.45 -1.24 -0.45 0.55 1.29]';
        ry=[0.30 1.15 1.46 1.09 0.20 -0.77 -1.39 -1.36 -0.69]';
        rz=zeros(length(rx),1);
        r=[rx,ry,rz];
%         figure;
%         plot(rx,ry,'o-');
%         hold on;
%         plot(rx_bar,ry_bar,'o-');
%         hold off;
        r= r*(lambda/2);
    else
        step=360/UCA;
        lambda=lightvel/Fc;
        thetas= (0:step:360-step);
        radius= (lambda/4)/ (sind(step/2));
        rx= radius*cosd(thetas);
        ry= radius*sind(thetas);
        rz=zeros(size(rx));
        r=[rx.',ry.',rz.'];
    end
           
end