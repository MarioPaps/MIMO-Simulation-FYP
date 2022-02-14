function [symbolsOut]=QPSKMod(bitsIn,radius,phi_rad)

    numsymbs=4;
    symbols= radius*exp(1i*(phi_rad+(0:numsymbs-1)*(pi/2)));
    symbolsOut=zeros(1,length(bitsIn)/2);
    pos00=[];
    pos01=[];
    pos11=[];
    pos10=[];
    for m=1:2:length(bitsIn)
       current=[bitsIn(m),bitsIn(m+1)];
       current= sprintf('%.0f',current);
       if(current=="00")
            pos00=[pos00 (m+1)/2];
       elseif(current=="01")
           pos01=[pos01 (m+1)/2];
       elseif(current=="11")
           pos11=[pos11 (m+1)/2];
       else
           pos10=[pos10 (m+1)/2];
       end
    end
    symbolsOut(pos00)=symbols(1);
    symbolsOut(pos01)=symbols(2);
    symbolsOut(pos11)=symbols(3);
    symbolsOut(pos10)=symbols(4);
end