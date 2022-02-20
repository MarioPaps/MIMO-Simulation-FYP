addpath('.\Utilities\');
addpath('.\Compute\');
%Define constants
NoSymbs=4000; %or 40
M=5; %# users
N_bar=16;
N=9;
Nsc=1;
Nc=31;
Ts=0.1e-6;
Tc=Ts*Nsc;
Tcs=Nc*Tc; 
Fc=20e9;
K=3; %3 paths per user
lightvel=3e8;
Fjvec=((0:Nsc-1)*(1/Tc))';
%% user data
%200 channel symbols for main user
bitstream=DataGen(2*NoSymbs);
radius=sqrt(2);
phi_rad= 43*(pi/180);
numsymbs=4;
symbols= radius*exp(1i*(phi_rad+(0:numsymbs-1)*(pi/2)));
A= QPSKMod(bitstream(1,:),radius,phi_rad); %desired user's symbols
MAI= zeros(M-1,NoSymbs);
MAIsyms=[0.1+1i;0.1-1i];
%Generate MAI
for user=1:M-1
    id1=randi([1,NoSymbs/2],1,NoSymbs/2);
    id2=randi([(NoSymbs/2)+1,NoSymbs],1,NoSymbs/2);
    
    MAI(user,id1)= MAIsyms(1);
    MAI(user,id2)= MAIsyms(2);
end
c=PNSeqGen();
ausers=[A,MAI(1,:),MAI(2,:),MAI(3,:),MAI(4,:)];
%% antenna array setup
[r,r_bar]=TxRxArr(lightvel,Fc);
%% Channel parameters
[delays,~,DODs,DOAs,VDops]= Channel_Param_Gen();
beta=[0.8,0.7,0.9; 0.02 0.04 0.15; 0.21 0.18 0.25; 0.35 0.31 0.3; 0.19 0.02 0.04];
beta=beta';
delays=round(delays/10);

f1j= computef(A,VDops(1,:),Fjvec,Fc,Tcs,lightvel,N_bar);
f2j= computef(MAI(1,:),VDops(2,:),Fjvec,Fc,Tcs,lightvel,N_bar);
f3j= computef(MAI(2,:),VDops(3,:),Fjvec,Fc,Tcs,lightvel,N_bar);
f4j= computef(MAI(3,:),VDops(4,:),Fjvec,Fc,Tcs,lightvel,N_bar);
f5j= computef(MAI(4,:),VDops(5,:),Fjvec,Fc,Tcs,lightvel,N_bar);

f=[f1j,f2j,f3j,f4j,f5j];
gamma1= computegamma(beta(:,1),DODs(1,:),Fjvec,r_bar);
gamma2= computegamma(beta(:,2),DODs(2,:),Fjvec,r_bar);
gamma3= computegamma(beta(:,3),DODs(3,:),Fjvec,r_bar);
gamma4= computegamma(beta(:,4),DODs(4,:),Fjvec,r_bar);
gamma5= computegamma(beta(:,5),DODs(5,:),Fjvec,r_bar);
gamma=[gamma1, gamma2,gamma3,gamma4,gamma5];

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
[Pn,lambda_min]= findPn(Rxx_prac,M);
del_est=[14,11,3]; 
vel_est=[20,66,120];
[cost1d]=OneDCost(del_est,vel_est,Fjvec,r,Pn,J,c(:,1),Nsc,K);
figure;
plot(20*log10(cost1d));
xlabel('DOA(degrees)'); ylabel('Gain(dB)');
[~,DOAest]=maxk(cost1d,K);





