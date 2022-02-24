addpath('.\Utilities\');
addpath('.\Compute\');
%Define constants
NoSymbs=1000;
M=2; %# users
N_bar=16;
N=9;
Nsc=5;
Nc=31;
Ts=0.1e-6;
Tc=Ts*Nsc;
Tcs=Nc*Tc; 
Fc=20e9;
K=3; %3 paths per user
lightvel=3e8;
Fjvec=((0:Nsc-1)*(1/Tc))';
%% Generate user and MAI data
radius=sqrt(2);
phi_rad= 43*(pi/180);
numsymbs=4;
symbols= radius*exp(1i*(phi_rad+(0:numsymbs-1)*(pi/2)));


phi_rad2= 0.05*(pi/180);
% bits=DataGen();
% filenames=["flamingo.jpg"];
% originalImage = imread(filenames(1));
% [rows,cols,~]= size(originalImage);
% P= rows*cols*24;
% for ind=1:length(filenames)
%     [bits,~,~]= fImageSource(filenames(ind),200);
% end
% bits=[bits;zeros(16,1)];
bitstream=DataGen(2*NoSymbs);
A= QPSKMod(bitstream(1,:),radius,phi_rad); %desired user's symbols
MAI= QPSKMod(bitstream,radius,phi_rad2);
Ptx= (1/length(A))* (A*A');
P2= (1/length(MAI))* (MAI*MAI');
a_i1= Demux(A(1,:),width(A),Nsc);
a_i2= Demux(MAI(1,:),width(A),Nsc);
disp('demux out');
%gold codes
c=PNSeqGen();
%% 
[r,r_bar]=TxRxArr(lightvel,Fc);

[delays,beta,DODs,DOAs,VDops]= Channel_Param_Gen();

f1j= computef(a_i1,VDops(1,:),Fjvec,Fc,Tcs,lightvel,N_bar);
f2j= computef(a_i2,VDops(2,:),Fjvec,Fc,Tcs,lightvel,N_bar);

f=[f1j,f2j];
gamma1= computegamma(beta(:,1:5),DODs(1,:),Fjvec,r_bar,K);
gamma2= computegamma(beta(:,6:10),DODs(2,:),Fjvec,r_bar,K);
gamma=[gamma1, gamma2];

ausers=[a_i1];
SNR_abs= 10^(20/10);
P_Tx=  (1/length(A))*( A*A');
P_MAI2= (1/width(MAI))* (MAI(1,:)* MAI(1,:)');
Pnoise= abs( P_Tx/SNR_abs);
%% H for equation 17
J = [zeros(1,2*Nc*Nsc-1) 0; eye(2*Nc*Nsc-1), zeros(2*Nc*Nsc-1,1)];
upper=NoSymbs/Nsc;
Next= 2*Nc*Nsc;
H=zeros(M*2*N*Nc*Nsc,K*Nsc);
for i=1:M
    for j=1:Nsc
        row_start= (i-1)*N*Next+1;
        row_end= i*N*Next;
        col_start= (j-1)*K+1;
        col_end= j*K;
        h1= DoppSTARmanifold(DOAs(i,1),delays(i,1),VDops(i,1),J,Fjvec(j),Nsc,r,c(:,i));
        h2= DoppSTARmanifold(DOAs(i,2),delays(i,2),VDops(i,2),J,Fjvec(j),Nsc,r,c(:,i));
        h3= DoppSTARmanifold(DOAs(i,3),delays(i,3),VDops(i,3),J,Fjvec(j),Nsc,r,c(:,i));
        H(row_start:row_end,col_start:col_end)=[h1,h2,h3];
    end
end
%% eqn17 compact
x=zeros(N*Next,upper);
for n=1:upper
    store= findX(ausers,f,gamma,H,J,1,Nsc,Nc,N,Next,K,n,upper,Pnoise);
    x(:,n)=store;
end
%% eqn 24
G= computeG(gamma,K);
Rxx_theor= covtheor(H,G,J,N,Nc,Nsc,M,Pnoise);
[Pn_theor,~]=findPn(Rxx_theor,M);
Rxx_prac= (1/width(x))* (x) * (x)';
predict=EstimateNumUsers(Rxx_prac,9,200)
%% 
h=[1 3; 2 4]
holder= h.*h
c= sum(h.*h,1)
Pn=ones(size(h));

