addpath('.\Utilities\');
addpath('.\Compute\');
%Define constants
NoSymbs=200; %or 40
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
%% 
c=PNSeqGen();
load("bitstream.mat");
bits=bitstream(1,:);
A= QPSKMod(bits,sqrt(2),deg2rad(43));
MAI= zeros(M-1,width(A));
bitsMAI= bitstream(2,:);
MAI(1,:)= QPSKMod(bitsMAI,1,deg2rad(0)); %symbol streams have the same power

for user=2:M-1
    MAI(user,:)= shuffle(MAI(1,:));
end
ai1=Demux(A,width(A),Nsc);
MAIres=cell(1,M-1);
for user=2:M
    %MAIres=[MAIres, Demux(MAI(user-1,:),width(A),Nsc)]; 
    MAIres{user-1}=Demux(MAI(user-1,:),width(A),Nsc);
end
ausers=[ai1, cell2mat(MAIres)];
ausers= unitymag(ausers); %every element with unity magnitude

%Tx outputs
[m1]=Tx(ai1,c(:,1),Nsc,N_bar,Tc);
[m2]=Tx(MAIres{1},c(:,2),Nsc,N_bar,Tc);
[m3]=Tx(MAIres{2},c(:,3),Nsc,N_bar,Tc);
[m4]=Tx(MAIres{3},c(:,4),Nsc,N_bar,Tc);
[m5]=Tx(MAIres{4},c(:,5),Nsc,N_bar,Tc);
mtot=[m1;m2;m3;m4;m5];
%testing what we actually transmit
atran=zeros(M*Nsc,M*width(A));
%% 
%Taking the average across chips of every symbol
for user=1:M
    rs= (user-1)*N_bar+1;
    re= user*N_bar;
    for n=1:NoSymbs
        cs=(n-1)*Nc+1;
        ce= n*Nc;
       check= (mtot(rs:re,cs:ce));
       atran(Nsc,(user-1)*n+n)=mean(mean(check));
    end

end

PTx= twodpower(m1);
PMAI= twodpower(m2);
%% 
% %% user data
% %200 channel symbols for main user
% bitstream=DataGen(2*NoSymbs);
% radius=sqrt(2);
% phi_rad= 43*(pi/180);
% numsymbs=4;
% symbols= radius*exp(1i*(phi_rad+(0:numsymbs-1)*(pi/2)));
% A= QPSKMod(bitstream(1,:),radius,phi_rad); %desired user's symbols
% MAI= zeros(M-1,NoSymbs);
% MAIsyms=[0.1+1i;0.1-1i];
% %Generate MAI
% for user=1:M-1
%     id1=randi([1,NoSymbs/2],1,NoSymbs/2);
%     id2=randi([(NoSymbs/2)+1,NoSymbs],1,NoSymbs/2);
%     
%     MAI(user,id1)= MAIsyms(1);
%     MAI(user,id2)= MAIsyms(2);
% end
% c=PNSeqGen();
% ausers=[A,MAI(1,:),MAI(2,:),MAI(3,:),MAI(4,:)];
%% Channel parameters
[r,r_bar]=TxRxArr(lightvel,Fc);
[delays,~,DODs,DOAs,VDops]= Channel_Param_Gen();
beta=[0.8*exp(1i*deg2rad(310)),0.7,0.9; 0.02 0.04 0.15; 0.21 0.18 0.25; 0.35 0.31 0.3; 0.19 0.02 0.04];
beta=beta.';
delays=round(delays/10);

f1j= computef(A,VDops(1,:),Fjvec,Fc,Tcs,lightvel,K);
f2j= computef(MAI(1,:),VDops(2,:),Fjvec,Fc,Tcs,lightvel,K);
f3j= computef(MAI(2,:),VDops(3,:),Fjvec,Fc,Tcs,lightvel,K);
f4j= computef(MAI(3,:),VDops(4,:),Fjvec,Fc,Tcs,lightvel,K);
f5j= computef(MAI(4,:),VDops(5,:),Fjvec,Fc,Tcs,lightvel,K);

f=[f1j,f2j,f3j,f4j,f5j];
gamma1= computegamma(beta(:,1),DODs(1,:),Fjvec,r_bar,K);
gamma2= computegamma(beta(:,2),DODs(2,:),Fjvec,r_bar,K);
gamma3= computegamma(beta(:,3),DODs(3,:),Fjvec,r_bar,K);
gamma4= computegamma(beta(:,4),DODs(4,:),Fjvec,r_bar,K);
gamma5= computegamma(beta(:,5),DODs(5,:),Fjvec,r_bar,K);
gamma=[gamma1, gamma2,gamma3,gamma4,gamma5];

SNR_abs= 10^(20/10);
PTx=twodpower(m1);
Pnoise= PTx/SNR_abs;
% P_Tx=  (1/length(A))*( A*A');
% P_MAI2= (1/width(MAI))* (MAI(1,:)* MAI(1,:)');
% P_MAI3= (1/width(MAI))* (MAI(2,:)* MAI(2,:)');
% P_MAI4= (1/width(MAI))* (MAI(3,:)* MAI(3,:)');
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
    store= findX(ausers,f,gamma,H,J,M,Nsc,Nc,N,Next,K,n,upper,0.1);
   %store= findX(atran,f,gamma,H,J,M,Nsc,Nc,N,Next,K,n,upper,0.1);
    x(:,n)=store;
end
% Pnoise= twodpower(x)/SNR_abs;
noise= sqrt(Pnoise/2)* (randn(size(x))+1i*randn(size(x)));
x=x+noise;
%% eqn 24
G= computeG(gamma,K);
Rxx_theor= covtheor(H,G,J,N,Nc,Nsc,M,Pnoise);
[Pn_theor,~]=findPn(Rxx_theor,M);
Rxx_prac= (1/width(x))* (x) * (x)';
%% 2d cost function inputs
[akj,Fkj]=findvecs(Fjvec,c(:,1),Nc,Nsc,Ts);
x_res= reshape(x,2*Nc*Nsc,[]);
Rxx_res=(1/width(x_res))* (x_res)*ctranspose(x_res);
[Pn_res,~]= findPn(Rxx_res,length(Rxx_res)-M);
%% 2d cost function
[cost2d,del_est,uk_est]= TwoDcost(N,Nc,Nsc,Fjvec,akj,Fkj,Pn_res,J);
figure;
rv_range=(1:140);
delay_range=(0:Nc*Nsc-1);
surf(rv_range,delay_range,20*log10((cost2d)),'FaceAlpha',1,'EdgeAlpha',0.5);
xlabel('Velocity(m/s)'); ylabel('Delay(Ts s)'); zlabel('Gain(dB)'); 
title('Joint Delay-Doppler Velocity Estimation');
%% 1d cost function
[Pn,lambda_min]= findPn(Rxx_theor,M);
del_est=[14,11,3]; 
vel_est=[20,66,120];
[cost1d]=OneDCost(del_est,vel_est,Fjvec,r,Pn_theor,J,c(:,1),Nsc,K);
[~,DOAest]=maxk(max(cost1d),K);
Colours = {'red','g','blue'};
figure;
for k=1:K
    txt = ['Multipath ',num2str(k) ];
    plot(20*log10(cost1d(k,:)),'Color',Colours{k},'DisplayName',txt);
    hold on;
end
xlabel('DOA(degrees)'); ylabel('Gain(dB)'); title('DOA Estimation'); legend('show');
hold off;
 %% gammakj amp est
k=1;
lambda_min= min(eig(Rxx_prac));
% [Pcompkjun,hkallj]=  gAmpSearch2(gamma,Rxx_prac,lambda_min,H,k,K,M,Nsc,N,Next,J,Fjvec,r,c(:,1));

%% try again
[Pcompkjun,hkallj]=  gAmpSearchN(gamma,Rxx_prac,lambda_min,H,k,K,M,Nsc,N,Next,J,Fjvec,r,c(:,1));

%% gamma_kj phase estimation for one path
phase_theor= angle(gamma(1,1:Nsc));
phi_rad= deg2rad(43);
phase_est=gPhSearch(x,f,ausers,Pcompkjun,hkallj,phi_rad,N,Next,Nsc,K);


 




