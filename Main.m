addpath('.\Utilities\');
addpath('.\Compute\');
%Define constants
NoSymbs=200;
M=5; %# users
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
c=PNSeqGen();
%% Generate user and MAI data
bits= round(rand(1,2*NoSymbs)) ; %bitstream(1,:); %;
A= QPSKMod(bits,sqrt(2),deg2rad(43));
ai1=Demux(A,width(A),Nsc);

MAI=zeros(M-1,length(A));
MAIres=cell(1,M-1);
rng default
for user=2:M
    bitsMAI= round(rand(1,2*NoSymbs));
    MAI(user-1,:)= QPSKMod(bitsMAI,0.2,deg2rad(0));
    MAIres{user-1}=Demux(MAI(user-1,:),width(A),Nsc);
end

ausers=[ai1, cell2mat(MAIres)];
ausers= unitymag(ausers); %every element with unity magnitude
ausers=ausers- real(mean(mean(ausers)));

%the symbol stream of user 1 has bigger power than the others
%Tx outputs
% [m1]=Tx(ai1,c(:,1),Nsc,N_bar,Tc);
% [m2]=Tx(MAIres{1},c(:,2),Nsc,N_bar,Tc);
% [m3]=Tx(MAIres{2},c(:,3),Nsc,N_bar,Tc);
% [m4]=Tx(MAIres{3},c(:,4),Nsc,N_bar,Tc);
% [m5]=Tx(MAIres{4},c(:,5),Nsc,N_bar,Tc);
% mtot=[m1;m2;m3;m4;m5];
% clear m1 m2 m3 m4 m5 
%% Define channel parameters
[r,r_bar]=TxRxArr(lightvel,Fc,"default");
N=length(r);
[delays,beta,DODs,DOAs,VDops]= Channel_Param_Gen(0,0);
%find f and gamma
f1j= computef(NoSymbs/Nsc,VDops(1,:),Fjvec,Fc,Tcs,lightvel,K);
f2j= computef(NoSymbs/Nsc,VDops(2,:),Fjvec,Fc,Tcs,lightvel,K);
f3j= computef(NoSymbs/Nsc,VDops(3,:),Fjvec,Fc,Tcs,lightvel,K);
f4j= computef(NoSymbs/Nsc,VDops(4,:),Fjvec,Fc,Tcs,lightvel,K);
f5j= computef(NoSymbs/Nsc,VDops(5,:),Fjvec,Fc,Tcs,lightvel,K);

f=[f1j,f2j,f3j,f4j,f5j]; clear f1j f2j f3j f4j f5j;
gamma1= computegamma(beta(:,1:5),DODs(1,:),Fjvec,r_bar,K);
gamma2= computegamma(beta(:,6:10),DODs(2,:),Fjvec,r_bar,K);
gamma3= computegamma(beta(:,11:15),DODs(3,:),Fjvec,r_bar,K);
gamma4= computegamma(beta(:,16:20),DODs(4,:),Fjvec,r_bar,K);
gamma5= computegamma(beta(:,21:25),DODs(5,:),Fjvec,r_bar,K);
gamma=[gamma1, gamma2,gamma3,gamma4,gamma5]; clear gamma1 gamma2 gamma3 gamma4 gamma5; 
G= computeG(gamma,K); %path power gain matrix

SNR_abs= 10^(20/10);
% P_Tx=  (1/length(A))*( A*A');
% %P_MAI2= (1/width(MAIsymbols))* (MAIsymbols(1,:)* MAIsymbols(1,:)');
% P_MAI2= (1/width(MAI))* (MAI(1,:)* MAI(1,:)');
% P_MAI3= (1/width(MAI))* (MAI(2,:)* MAI(2,:)');
% P_MAI4= (1/width(MAI))* (MAI(3,:)* MAI(3,:)');
% P_MAI5= (1/width(MAI))* (MAI(4,:)* MAI(4,:)');
Pnoise= abs(1/SNR_abs);
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
%% find x-equation 17
x=zeros(N*Next,L);
% x=[];
tic;
for n=1:L
    x(:,n)=findX(ausers,f,gamma,H,J,M,Nsc,Nc,N,Next,K,n,L,Pnoise);
end
% load("x.mat");
%noise= sqrt(Pnoise/2)* (randn(size(x))+1i*randn(size(x)));
%x=x+noise;
toc;
%% eqn 24
%load("x.mat");

Rxx_theor= covtheor(H,G,J,N,Nc,Nsc,M,Pnoise);
[Pn_theor,~]=findPn(Rxx_theor,M*Nsc);
Rxx_prac= (1/L)* (x) * (x)';
%% Delay-Velocity cost function inputs
[akj,Fkj]=findvecs(Fjvec,(1:140),c(:,1),Nc,Nsc,Ts);
x_res= reshape(x,2*Nc*Nsc,[]);
Rxx_res= (1/width(x_res))* (x_res)*ctranspose(x_res);
[Pn_res,~]= findPn(Rxx_res,length(Rxx_res)-M*Nsc);
%% Delay-Velocity cost function 
tic;
[cost2d,del_est,uk_est]= faster2dcost(K,Nc,Nsc,delays(1,:),(1:140),akj,Fkj,Pn_res,J);
%[cost2d,del_est,uk_est]= TwoDcost(K,Nc,Nsc,delays(1,:),(1:140),akj,Fkj,Pn_res,J);
figure;
surf((1:140),(0:Nc*Nsc-1),10*log((cost2d)),'FaceAlpha',1,'EdgeAlpha',0.5);
xlabel('Velocity(m/s)'); ylabel('Delay(Ts s)'); zlabel('Gain(dB)');
title('Joint Delay-Doppler Velocity Estimation');
shading('interp');
colormap('jet');
toc;
%% DOA Cost Function
[Pn,~]= findPn(Rxx_theor,M*Nsc);
del_est=[140,110,30]; 
vel_est=[20,66,120];
%[cost1d]=OneDCost(del_est,uk_est,Fjvec,r,Pn,J,c(:,1),Nsc,K);
[cost1d]=experiment(del_est,vel_est, (1:360), Fjvec,r,Pn,J,c(:,1),Nsc,K);
Colours = {'red','g','blue'};
figure;
for k=1:K
    txt = ['Multipath ',num2str(k) ];
    plot(20*log10(cost1d(k,:)),'Color',Colours{k},'DisplayName',txt);
    hold on;
end
xlabel('DOA(degrees)'); ylabel('Gain(dB)'); title('DOA Estimation'); legend('show');
DOAest= findMaxofPath(cost1d);
hold off;

%% gamma_kj ampl. estimation for one path
k=1;
lambda_min= min(eig(Rxx_prac));
%[Pcompkjun,hkallj]= gAmpSearch2(gamma,Rxx_prac,lambda_min,H,k,K,M,Nsc,N,Next);
[Pcompkjun,hkallj]=  gAmpSearch2(gamma,Rxx_prac,lambda_min,H,k,K,M,Nsc,N,Next,J,Fjvec,r,c(:,1));
%% try again
[Pcompkjun,hkallj]=  gAmpSearchN(gamma,Rxx_prac,lambda_min,H,k,K,M,Nsc,N,Next,J,Fjvec,r,c(:,1));

%% gamma_kj phase estimation for one path
phase_theor= angle(gamma(1,1:Nsc));
phi_rad= deg2rad(43);
phase_est=gPhSearch(x,f,ausers,Pcompkjun,hkallj,phi_rad,N,Next,Nsc,K);


% gammakj_desuser= (atan2(gamma(1,1:Nsc))); %wrapTo2Pi
% display_res= [];
% for iter=1:length(phase_est)
%     display_res=[display_res; gammakj_desuser(iter),phase_est(iter)];
% end
% figure;
% bar(display_res);
% xlabel('Subcarrier Index'); ylabel('Phase Estimate(rad)');

%% faster 2dcost
% tic;
% [cost2d,del_est,uk_est]= faster2dcost(K,Nc,Nsc,delays(1,:),akj,Fkj,Pn_res,J);
% figure;
% rv_range=(1:140);
% delay_range=(0:Nc*Nsc-1);
% surf(rv_range,delay_range,20*log10((cost2d))-40,'FaceAlpha',1,'EdgeAlpha',0.5);
% xlabel('Velocity(m/s)'); ylabel('Delay(Ts s)'); zlabel('Gain(dB)');
% title('Joint Delay-Doppler Velocity Estimation');
% toc;
%% playing with tall arrays
y=dot(x(:,1),x(:,1));





