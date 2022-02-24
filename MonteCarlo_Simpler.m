%constants
addpath('.\Utilities\');
addpath('.\Compute\');
%Define constants
NoSymbs=2000; 
M=5; %# users
N_bar=16;
N=9;
Nsc=1;
Nc=31;
Ts=0.1e-6;
Tc=Ts*Nsc;
Tcs=Nc*Tc; 
Fc=20e9;
K=3; %2 paths per user
lightvel=3e8;
Fjvec=((0:Nsc-1)*(1/Tc))';
%% inputs
len=2*NoSymbs*Nsc;
bits= DataGen(len);
radius=sqrt(2);
phi_rad= 43*(pi/180);
numsymbs=4;
symbols= radius*exp(1i*(phi_rad+(0:numsymbs-1)*(pi/2)));
A= QPSKMod(bits,radius,phi_rad); %desired user's symbols
MAI= zeros(M-1,width(A));
MAIsyms=[0.1+1i;0.1-1i];
%Generate MAI
for user=1:M-1
    id1=randi([1,width(A)/2],1,width(A)/2);
    id2=randi([(width(A)/2)+1,width(A)],1,width(A)/2);
    MAI(user,id1)= MAIsyms(1);
    MAI(user,id2)= MAIsyms(2);
end
a_i1= Demux(A(1,:),width(A),Nsc);
a_i2= Demux(MAI(1,:),width(A),Nsc);
a_i3= Demux(MAI(2,:),width(A),Nsc);
a_i4= Demux(MAI(3,:),width(A),Nsc);
a_i5= Demux(MAI(3,:),width(A),Nsc);
disp('demux out');
%gold codes
c=PNSeqGen();
%% channel
[r,r_bar]=TxRxArr(lightvel,Fc);
[delays,~,DODs,DOAs,VDops]= Channel_Param_Gen();
beta=[0.8*exp(1i*deg2rad(310)),0.7,0.9; 0.02 0.04 0.15; 0.21 0.18 0.25; 0.35 0.31 0.3; 0.19 0.02 0.04];
beta=beta';
delays=round(delays/10);

f1j= computef(A,VDops(1,:),Fjvec,Fc,Tcs,lightvel,N_bar);
f2j= computef(MAI(1,:),VDops(2,:),Fjvec,Fc,Tcs,lightvel,N_bar);
f3j= computef(MAI(2,:),VDops(3,:),Fjvec,Fc,Tcs,lightvel,N_bar);
f4j= computef(MAI(3,:),VDops(4,:),Fjvec,Fc,Tcs,lightvel,N_bar);
f5j= computef(MAI(4,:),VDops(5,:),Fjvec,Fc,Tcs,lightvel,N_bar);

f=[f1j,f2j,f3j,f4j,f5j];
gamma1= computegamma(beta(:,1),DODs(1,:),Fjvec,r_bar,K);
gamma2= computegamma(beta(:,2),DODs(2,:),Fjvec,r_bar,K);
gamma3= computegamma(beta(:,3),DODs(3,:),Fjvec,r_bar,K);
gamma4= computegamma(beta(:,4),DODs(4,:),Fjvec,r_bar,K);
gamma5= computegamma(beta(:,5),DODs(5,:),Fjvec,r_bar,K);
gamma=[gamma1, gamma2,gamma3,gamma4,gamma5];
ausers=[a_i1,a_i2,a_i3,a_i4,a_i5];

SNR_abs= 10^(20/10);
P_Tx=  (1/length(A))*( A*A');
P_MAI2= (1/width(MAI))* (MAI(1,:)* MAI(1,:)');
P_MAI3= (1/width(MAI))* (MAI(2,:)* MAI(2,:)');
P_MAI4= (1/width(MAI))* (MAI(3,:)* MAI(3,:)');

%% H for equation 17
J = [zeros(1,2*Nc*Nsc-1) 0; eye(2*Nc*Nsc-1), zeros(2*Nc*Nsc-1,1)];
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
%% 
stepsize=0.1;
thetas=(1:stepsize:360);
PTx= (1/length(A))*(A*A');
L=200; %200 snapshots
SNR=5:5:30; %5-30 dB
SNR=10.^(SNR/10);
numtrials=100;
upper=width(a_i1);
RMSEradvel=zeros(numtrials,length(SNR));
RMSEgamma= zeros(numtrials,length(SNR));
RMSEDOA= zeros(numtrials,length(SNR));
%% obtain noiseless x
x=zeros(N*Next,upper);
for n=1:upper
       x(:,n)=findX(ausers,f,gamma,H,J,M,Nsc,Nc,N,Next,K,n,upper,0.1);
end
%% Trial 1 - looping through SNRs
[akj,Fkj]=findvecs(Fjvec,c(:,1),Nc,Nsc,Ts);
for ind=1:length(SNR)
    Pnoise=PTx/SNR(ind);
    noise= sqrt(Pnoise/2)* (randn(size(x))+1i*randn(size(x)));
    x_noisy= x+noise;
    x_res= reshape(x_noisy,2*Nc*Nsc,[]);
    Rxx_res=(1/width(x_res))* (x_res)*ctranspose(x_res);
    Rxx_prac= (1/width(x_noisy))* (x_noisy) * (x_noisy)';
    [Pn,lambda_min]= findPn(Rxx_prac,M);
    [Pn_res,~]= findPn(Rxx_res,M);

    %rmse of radial velocity
    [cost2d,del_est,uk_est]= TwoDcost(K,Nc,Nsc,delays(1,:),akj,Fkj,Pn_res,J);
%     uk_est(3)=118;
    RMSEradvel(1,ind)= findRMSE(VDops(1,:),uk_est);
%% 

% %     %rmse of doa
%     [cost1d]=experiment(del_est,uk_est,thetas,Fjvec,r,Pn,J,c(:,1),Nsc,K);
%     DOAest= findMaxofPath(cost1d);
%     DOAest= thetas(DOAest);
%     RMSEDOA(1,ind)= findRMSE(DOAs(1,:),DOAest);

    disp(ind);
end


