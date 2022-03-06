function[unity]=unitymag(in)
    [r,c]=size(in);
    unity=zeros(size(in));
    for ind=1:c
        for knd=1:r
            unity(knd,ind)= in(knd,ind)/ abs(in(knd,ind));
        end
    end
end