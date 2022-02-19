addpath('.\Utilities\');
addpath('.\Compute\');
%Define constants
NoSymbs=200;
M=4; %# users
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
%create 400 bits long bitstream -> 200 channel symbols via QPSK modulation
load("bitstream.mat");

radius=sqrt(2);
phi_rad= 43*(pi/180);
numsymbs=4;
symbols= radius*exp(1i*(phi_rad+(0:numsymbs-1)*(pi/2)));

MAI= zeros(M-1,NoSymbs);
MAIsyms=[0.1+1i;0.1-1i];
%Generate MAI
for user=1:M-1
    id1=randi([1,NoSymbs/2],1,NoSymbs/2);
    id2=randi([(NoSymbs/2)+1,NoSymbs],1,NoSymbs/2);
    
    MAI(user,id1)= MAIsyms(1);
    MAI(user,id2)= MAIsyms(2);
end
% bits=DataGen();
% filenames=["flamingo.jpg"];
% originalImage = imread(filenames(1));
% [rows,cols,~]= size(originalImage);
% P= rows*cols*24;
% for ind=1:length(filenames)
%     [bits,~,~]= fImageSource(filenames(ind),200);
% end
% bits=[bits;zeros(16,1)];
A= QPSKMod(bitstream(1,:),radius,phi_rad); %desired user's symbols

%demultiplexer output for every user
%The 200 symbols of every user are reshaped into a 5*col matrix
%the interfering users go through a different channel with BPSK modulation
a_i1= Demux(A(1,:),NoSymbs,Nsc);
a_i2= Demux(MAI(1,:),NoSymbs,Nsc);
a_i3= Demux(MAI(2,:),NoSymbs,Nsc);
a_i4= Demux(MAI(3,:),NoSymbs,Nsc);

disp('demux out');
%gold codes
c=PNSeqGen();
%% antenna array setup
[r,r_bar]=TxRxArr(lightvel,Fc);
%% Define channel parameters
[delays,beta,DODs,DOAs,VDops]= Channel_Param_Gen();

f1j= computef(a_i1,VDops(1,:),Fjvec,Fc,Tcs,lightvel,N_bar);
f2j= computef(a_i2,VDops(2,:),Fjvec,Fc,Tcs,lightvel,N_bar);
f3j= computef(a_i3,VDops(3,:),Fjvec,Fc,Tcs,lightvel,N_bar);
f4j= computef(a_i4,VDops(4,:),Fjvec,Fc,Tcs,lightvel,N_bar);
%f5j= computef(a_i5,VDops(5,:),Fjvec,Fc,Tcs,lightvel,N_bar);

f=[f1j,f2j,f3j,f4j];
gamma1= computegamma(beta(:,1:5),DODs(1,:),Fjvec,r_bar);
gamma2= computegamma(beta(:,6:10),DODs(2,:),Fjvec,r_bar);
gamma3= computegamma(beta(:,11:15),DODs(3,:),Fjvec,r_bar);
gamma4= computegamma(beta(:,16:20),DODs(4,:),Fjvec,r_bar);
%gamma5= computegamma(beta(:,21:25),DODs(5,:),Fjvec,r_bar);
gamma=[gamma1, gamma2,gamma3,gamma4];

ausers=[a_i1,a_i2,a_i3,a_i4];
SNR_abs= 10^(20/10);
P_Tx=  (1/length(A))*( A*A');
P_MAI2= (1/width(MAI))* (MAI(1,:)* MAI(1,:)');
P_MAI3= (1/width(MAI))* (MAI(2,:)* MAI(2,:)');
P_MAI4= (1/width(MAI))* (MAI(3,:)* MAI(3,:)');
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
    store= findX(ausers,f,gamma,H,J,M,Nsc,Nc,N,Next,K,n,upper,Pnoise);
    x(:,n)=store;
end
%% eqn 24
G= computeG(gamma,K);
Rxx_theor= covtheor(H,G,J,N,Nc,Nsc,M,Pnoise);
[Pn_theor,~]=findPn(Rxx_theor,1);
Rxx_prac= (1/width(x))* (x) * (x)';
%% 2d cost function inputs
[akj,Fkj]=findvecs(Fjvec,c(:,1),Nc,Nsc,Ts);
x_res= reshape(x,2*Nc*Nsc,[]);
Rxx_res=(1/width(x_res))* (x_res)*ctranspose(x_res);
[Pn_res,~]= findPn(Rxx_res,M);
%% 2d cost function
[cost2d,del_est,uk_est]= TwoDcost(N,Nc,Nsc,Fjvec,akj,Fkj,Pn_res,J);
figure;
surf(20*log10((cost2d)),'FaceAlpha',1,'EdgeAlpha',0.5);
xlabel('Velocity(m/s)'); ylabel('Delay(Ts s)'); zlabel('Gain(dB)');
%% 1d cost function
[Pn,lambda_min]= findPn(Rxx_theor,M);
del_est=[141,111,31]; %zero-based indexing
vel_est=[20,66,120];
[cost1d]=OneDCost(del_est,vel_est,Fjvec,r,Pn,J,c(:,1),Nsc,K);
figure;
plot(20*log10(cost1d));
xlabel('DOA(degrees)'); ylabel('Gain(dB)');
[~,DOAest]=maxk(cost1d,K);
%% gamma_kj ampl. estimation for one path
k=1;
[Pcompkjun,hkallj]= gAmpSearch(gamma,Rxx_prac,H,k,M,Nsc,N,Next);
psi0=phi_rad;
%% gamma_kj phase estimation for one path
phase_est=gPhSearch(x_noisy,f,ausers,Pcompkjun,hkallj,psi0,N,Next,Nsc);
gammakj_desuser= wrapTo2Pi(angle(gamma(1,1:Nsc)));
display_res= [];
for iter=1:length(phase_est)
    display_res=[display_res; gammakj_desuser(iter),phase_est(iter)];
end
figure;
bar(display_res);
xlabel('Subcarrier Index'); ylabel('Phase Estimate(rad)');
