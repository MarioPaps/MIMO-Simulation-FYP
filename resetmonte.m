addpath('.\Utilities\');
addpath('.\Compute\');
%Define constants
NoSymbs=1000;
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
radius=sqrt(2);
phi_rad= 43*(pi/180);
bits=DataGen(2*NoSymbs);
A= QPSKMod(bits,radius,phi_rad);
bitsMAI=DataGen(2*NoSymbs);
MAI= zeros(M-1,width(A));
MAI(1,:)= QPSKMod(bitsMAI,radius,phi_rad);

for user=2:M-1
    MAI(user,:)=shuffle(MAI(1,:));
end
a_i1= Demux(A(1,:),width(A),Nsc);
a_i2= Demux(MAI(1,:),width(A),Nsc);
a_i3= Demux(MAI(2,:),width(A),Nsc);
a_i4= Demux(MAI(3,:),width(A),Nsc);
disp('demux out');
%gold codes
c=PNSeqGen();
%% Define channel parameters
[r,r_bar]=TxRxArr(lightvel,Fc);
[delays,beta,DODs,DOAs,VDops]= Channel_Param_Gen();

f1j= computef(a_i1,VDops(1,:),Fjvec,Fc,Tcs,lightvel,K);
f2j= computef(a_i2,VDops(2,:),Fjvec,Fc,Tcs,lightvel,K);
f3j= computef(a_i3,VDops(3,:),Fjvec,Fc,Tcs,lightvel,K);
f4j= computef(a_i4,VDops(4,:),Fjvec,Fc,Tcs,lightvel,K);


f=[f1j,f2j,f3j,f4j];
gamma1= computegamma(beta(:,1:5),DODs(1,:),Fjvec,r_bar,K);
gamma2= computegamma(beta(:,6:10),DODs(2,:),Fjvec,r_bar,K);
gamma3= computegamma(beta(:,11:15),DODs(3,:),Fjvec,r_bar,K);
gamma4= computegamma(beta(:,16:20),DODs(4,:),Fjvec,r_bar,K);
%gamma5= computegamma(beta(:,21:25),DODs(5,:),Fjvec,r_bar,K);
gamma=[gamma1, gamma2,gamma3,gamma4];

ausers=[a_i1,a_i2,a_i3,a_i4];
ausers= unitymag(ausers);
P_Tx=  (1/length(A))*( A*A');
P_MAI2= (1/width(MAI))* (MAI(1,:)* MAI(1,:)');
P_MAI3= (1/width(MAI))* (MAI(2,:)* MAI(2,:)');
P_MAI4= (1/width(MAI))* (MAI(3,:)* MAI(3,:)');
%% H for equation 17
J = [zeros(1,2*Nc*Nsc-1) 0; eye(2*Nc*Nsc-1), zeros(2*Nc*Nsc-1,1)];
L=NoSymbs/Nsc;
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
SNR_db=5:5:30; %5-30 dB
SNR=10.^(SNR_db/10);
numtrials=100;
RMSEradvel=zeros(numtrials,length(SNR));
RMSEgamma= zeros(numtrials,length(SNR));
RMSEDOA= zeros(numtrials,length(SNR));
%% obtain noiseless x
x=zeros(N*Next,L);
for n=1:L
       x(:,n)=findX(ausers,f,gamma,H,J,M,Nsc,Nc,N,Next,K,n,L,0.1);
end
%% obtain x for every trial and SNR 
x_noise= cell(numtrials,length(SNR));
for trial=1:numtrials
    for ind=1:length(SNR)
        Pnoise=P_Tx/SNR(ind);
        noise= sqrt(Pnoise/2)* (randn(size(x))+1i*randn(size(x)));
        x_noise{trial,ind}=x+noise;
    end
end
 %% obtain uk estimates for every SNR

 [akj,Fkj]=findvecs(Fjvec,c(:,1),Nc,Nsc,Ts);
 vel_est=[];
tic;
 for trial=2:2
      uk_est=zeros(length(SNR),K);
     for ind=1:length(SNR)
        x_res= reshape(x_noise{trial,ind},2*Nc*Nsc,[]);
        Rxx_res= (1/width(x_res))* (x_res)*ctranspose(x_res);
        [Pn_res,~]= findPn(Rxx_res,length(Rxx_res)-M);
        [cost2d,del_est,uk_est(ind,:)]= TwoDcost(K,Nc,Nsc,delays(1,:),akj,Fkj,Pn_res,J);
        RMSEradvel(trial,ind)=findRMSE(VDops(1,:),uk_est(ind,:));
        
     end
     vel_est=[vel_est,uk_est];
 end
 toc;
 %% obtain DOA estimates for every SNR
 tic;
 for trial=1:1
     DOA_est=zeros(length(SNR),K);
    for ind=1:length(SNR)
        uk_est= vel_est(:, (trial-1)*K+1: trial*K);
        xcurr= x_noise{trial,ind};
        Rxx_prac= (1/L)* (xcurr) * ctranspose(xcurr);
        [Pn,lambda_min]= findPn(Rxx_prac,M);
        [cost1d]=experiment(del_est,uk_est,thetas,Fjvec,r,Pn,J,c(:,1),Nsc,K);
        DOAest= findMaxofPath(cost1d);
        DOAest= thetas(DOAest);
        RMSEDOA(trial,ind)= findRMSE(DOAs(1,:),DOAest);
    end
 end
toc;
 %% Plot results
xaxis= SNR_db*L;
figure;
stem(mean(RMSEradvel));
 
%    [cost2d,del_est,uk_est]= TwoDcost(K,Nc,Nsc,delays(1,:),akj,Fkj,Pn_res,J);
% %     uk_est(3)=118;
%     RMSEradvel(1,ind)= findRMSE(VDops(1,:),uk_est);

