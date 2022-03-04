function [symbolsOut]=BPSKmod(bitsIn)

    numsymbs=2;
    symbols= [1,-1];
    symbolsOut=zeros(1,length(bitsIn)/2);
    pos0=[];
    pos1=[];
    for m=1:length(bitsIn)
       current=bitsIn(m);
       current= sprintf('%.0f',current);
       if(current=="0")
            pos0=[pos0 m];
       else
            pos1=[pos1 m];
       end
    end
    symbolsOut(pos0)=symbols(1);
    symbolsOut(pos1)=symbols(2);
end