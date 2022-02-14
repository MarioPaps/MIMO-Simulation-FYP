function[power]= twodpower(in)
    [Nbar,len]=size(in);
    acc=0;
    for ind=1:len
        singlecolpow= (1/Nbar)*ctranspose(in(:,ind))*in(:,ind);
        acc=acc+singlecolpow;
    end

    power=acc;
end