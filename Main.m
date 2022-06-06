addpath('.\Utilities\');
addpath('.\Compute\');
%Define constants
L=40; %number of snaphots of x[n]
Nsc=5;
NoSymbs= L*Nsc;
M=5; %# users
N_bar=16;
Nsc=1;
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
bits= round(rand(1,2*NoSymbs)) ; %
% load("bitstream.mat");
% bits=bitstream(1,:);
A= QPSKMod(bits,sqrt(2),deg2rad(43));
ai1=Demux(A,width(A),Nsc);

MAI=zeros(M-1,length(A));
MAIres=cell(1,M-1);

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
[r,r_bar]=TxRxArr(lightvel,Fc,9);
N=length(r);
[delays,beta,DODs,DOAs,VDops]= ChannelParam(0,0,0,Nsc);
delays(1,:)=[42,85,103];
VDops(1,:)=[25,72,106];
DOAs(1,:)=[62,150,220];
%beta=real(beta);

%find f and gamma
f1j= computef(NoSymbs/Nsc,VDops(1,1:K),Fjvec,Fc,Tcs,lightvel,K);
f2j= computef(NoSymbs/Nsc,VDops(2,1:K),Fjvec,Fc,Tcs,lightvel,K);
f3j= computef(NoSymbs/Nsc,VDops(3,1:K),Fjvec,Fc,Tcs,lightvel,K);
f4j= computef(NoSymbs/Nsc,VDops(4,1:K),Fjvec,Fc,Tcs,lightvel,K);
f5j= computef(NoSymbs/Nsc,VDops(5,1:K),Fjvec,Fc,Tcs,lightvel,K);

f=[f1j,f2j,f3j,f4j,f5j]; clear f1j f2j f3j f4j f5j;
gamma1= computegamma(beta(:,1:Nsc),DODs(1,:),Fjvec,r_bar,K);
gamma2= computegamma(beta(:,Nsc+1:2*Nsc),DODs(2,:),Fjvec,r_bar,K);
gamma3= computegamma(beta(:,2*Nsc+1:3*Nsc),DODs(3,:),Fjvec,r_bar,K);
gamma4= computegamma(beta(:,3*Nsc+1:4*Nsc),DODs(4,:),Fjvec,r_bar,K);
gamma5= computegamma(beta(:,4*Nsc+1:5*Nsc),DODs(5,:),Fjvec,r_bar,K);
gamma=[gamma1, gamma2,gamma3,gamma4,gamma5];
gamma=real(gamma);
G= computeG(gamma,K);
clear gamma1 gamma2 gamma3 gamma4 gamma5; 

SNR_abs= 10^(20/10);
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
%         for k=1:K
%             h= DoppSTARmanifold(DOAs(i,k),delays(i,k),VDops(i,k),J,Fjvec(j),Nsc,r,c(:,i));
%         end
        h1= DoppSTARmanifold(DOAs(i,1),delays(i,1),VDops(i,1),J,Fjvec(j),Nsc,r,c(:,i));
        h2= DoppSTARmanifold(DOAs(i,2),delays(i,2),VDops(i,2),J,Fjvec(j),Nsc,r,c(:,i));
        h3= DoppSTARmanifold(DOAs(i,3),delays(i,3),VDops(i,3),J,Fjvec(j),Nsc,r,c(:,i));
        H(row_start:row_end,col_start:col_end)=[h1,h2,h3];
      % H(row_start:row_end,col_start:col_end)=h1;
    end
end
clear h1 h2 h3;
%% find x-equation 17
x=zeros(N*Next,L);
% x=[];
tic;
for n=1:L
    x(:,n)=findX(ausers,f,gamma,H,J,M,Nsc,Nc,N,Next,K,n,L,Pnoise);
end
% load("x.mat");
noise= sqrt(Pnoise/2)* (randn(size(x))+1i*randn(size(x)));
x=x+noise;
toc;
%% Delay-Velocity cost function inputs
Rxx_prac= (1/L)* (x) * (x)';
[akj,Fkj]=findvecs(Fjvec,(1:140),c(:,1),Nc,Nsc,Ts);
x_res= reshape(x,2*Nc*Nsc,[]);
Rxx_res= (1/width(x_res))* (x_res)*ctranspose(x_res);
[Pn_res,~]= findPn(Rxx_res,length(Rxx_res)-M*Nsc);
%% Delay-Velocity cost function 
tic;
[cost2d,del_est,uk_est]= faster2dcost(K,Nc,Nsc,delays(1,:),(1:140),akj,Fkj,Pn_res,J);
figure;
surf((1:140),(0:Nc*Nsc-1),15*log10((cost2d))-10,'FaceAlpha',1,'EdgeAlpha',0.5);
xlabel('Velocity(m/s)'); ylabel('Delay(Ts s)'); zlabel('Equation NUM Amplitude(dB)');
zlim([0 150]);
title('Joint Delay-Doppler Velocity Estimation');
shading('interp');
colormap('jet');
ax = gca; 
ax.FontSize = 11; 
toc;
%% DOA Cost Function
Rxx_theor= covtheor(H,G,N,Nc,K,Nsc,M,Pnoise);
[Pn,~]= findPn(Rxx_theor,M*Nsc);
del_est=delays(1,:); 
vel_est=VDops(1,:);
%[cost1d]=OneDCost(del_est,uk_est,Fjvec,r,Pn,J,c(:,1),Nsc,K);
[cost1d]=experiment(del_est,vel_est, (1:360), Fjvec,r,Pn,J,c(:,1),Nsc,K);
Colours = {'red','g','blue'};
figure;
for k=1:K
    txt = ['Multipath ',num2str(k) ];
    plot(20*log10(cost1d(k,:)),'Color',Colours{k},'DisplayName',txt);
    hold on;
end
grid on;
xlabel('DOA(degrees)'); ylabel('Equation NUM Amplitude(dB)'); title('DOA Estimation'); legend('show');
DOAest= findMaxofPath(cost1d);
hold off;
%% gamma_kj ampl. estimation for all paths on 1 subcarrier
gammapred=zeros(1,K);
Hj=H(1:height(H)/5,:);
for k=1:K
    gammapred(k)=gammaAmpSearch(Rxx_prac,gamma,Hj,k,N,Nc,Nsc);
end
gammapred= fix(gammapred.*100)./100;
gamma_actual= fix(gamma(1:K,1).*100)./100;
%% bar char plotting
plotmat=zeros(K,2);
for k=1:K
    plotmat(k,:)=[gammapred(k),gamma_actual(k)];
end
figure;
h=bar(plotmat);
xlabel('Multipath Index');ylabel('Fading Coefficient Amplitude');
title('Fading Coefficient Amplitude Estimation');
ax = gca; 
ax.FontSize = 11; 

%gamma_est=gammaAmpSearch(Rxx_prac,gamma,Hj,k,N,Nc,Nsc);
%gamma_est= fix(gamma_est*100)/100
%% 

%gamma_est=gampsearchN(Rxx_prac,lambda_min,Hj,gamma,Gj,N,Nc,Nsc,Next);
%[Pcompkjun,hkallj]=  gAmpSearchN(gamma,Rxx_prac,lambda_min,H,k,K,M,Nsc,N,Next,J,Fjvec,r,c(:,1));
% phase_est=gPhSearch(x,f,ausers,Pcompkjun,hkallj,deg2rad(43),N,Next,Nsc,K);
%% 
dd = [1.5818 1.5332; 1 1];
bar(dd)

%% 






