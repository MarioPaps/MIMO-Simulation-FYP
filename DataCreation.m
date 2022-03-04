addpath('.\Utilities\');
addpath('.\Compute\');
filenames=["flamingo.jpg"];

%define constants

radius=sqrt(2);
phi_rad= 43*(pi/180);

P= 4*4*24;

[bitsOut,x,y]= fImageSource(filenames(1),400);
A= QPSKMod(bitsOut,sqrt(2),phi_rad);
%% 

%other users
phi=0.663225; %38 degrees = 0.663225 radians
id=fopen('input1.txt','r');
formatSpec='%c';
A2=fscanf(id,formatSpec);
ascivals=double(A2);
[~ , cols]=size(ascivals);
encoder=strings(1,cols);
s1='0';
for i=1:1:cols
    temp=ascivals(i);
    temp=dec2bin(temp);
    if(length(temp)~=8)
        for j=1:(8-length(temp))
            temp=append(s1,temp);
        end
    end
    encoder(1,i)=temp;
end
%% 
MAI= encoder(1:25);
temp=[];
MAIbits=[];
for ind=1:length(MAI)
    temp=MAI(ind);
    sol = split(temp,"")
    keep= sol(2:9)
    keep= str2double(keep)
    keep= keep';
    MAIbits=[MAIbits,keep];
end
MAIsymbols= BPSKmod(MAIbits);
