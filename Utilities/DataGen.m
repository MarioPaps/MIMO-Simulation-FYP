function[bitstream]= DataGen()
    len=4000;
    num_ones=len/2;
    bitstream= zeros(1,len);

 
   indices=randperm(len,num_ones);
   temp_stream(indices)=1;
   bitstream(1,:)= temp_stream;

end