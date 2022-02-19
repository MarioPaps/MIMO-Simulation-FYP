function[bitstream]= DataGen(len)
    num_ones=len/2;
    bitstream= zeros(1,len);

   indices=randperm(len);
   %assign first half of indices to 1
   pos_one= indices(1:(length(indices)/2));
   %assign second half of indices to 0
   pos_zero=indices((length(indices)/2)+1: length(indices));
   bitstream(pos_one)=1;
   bitstream(pos_zero)=0;
end