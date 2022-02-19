%constants
NoSymbs=200;
Fc=20e9;
Ts=0.1e-6;
Nsc=5;
K=2;
N_bar=16;
N=9;
M=4;
%% 4 users - equipowered
len=NoSymbs*Nsc;
bits= DataGen(len);
radius=sqrt(2);
phi_rad= 43*(pi/180);
numsymbs=4;
symbols= radius*exp(1i*(phi_rad+(0:numsymbs-1)*(pi/2)));
A= QPSKMod(bits,radius,phi_rad); %desired user's symbols
MAI= zeros(M-1,len);
MAIsyms=[0.1+1i;0.1-1i];
%Generate MAI
for user=1:M-1
    id1=randi([1,NoSymbs/2],1,NoSymbs/2);
    id2=randi([(NoSymbs/2)+1,NoSymbs],1,NoSymbs/2);
    MAI(user,id1)= MAIsyms(1);
    MAI(user,id2)= MAIsyms(2);
end
%% 
