%Objective: Generate  symbols per user that will be sent to the demux prior
%to the transmitter
M=5; %# users
%user 1 is the desired user
K=3; %3 paths per user

%create 400 bits long bitstream -> 200 symbols going into demux
len=400;
num_ones=len/2;
bitstream= zeros(M,len);
temp_stream=zeros(1,len);

for user=1:M
   indices=randperm(len,num_ones);
   temp_stream(indices)=1;
   bitstream(user,:)= temp_stream;
end

save('bitstream.mat','bitstream');