addpath('.\Utilities\');
addpath('.\Compute\');
addpath('.\Beampattern_files\');
M=5; %# users   
N_bar=16;
L=40; %number of snaphots of x[n]
Nsc=1;
NoSymbs= L*Nsc;
Nc=31;
Ts=0.1e-6;
Tc=Ts*Nsc;
Tcs=Nc*Tc; 
Fc=20e9;
K=3; %3 paths per user
lightvel=3e8;
lambda= lightvel/Fc;
Fjvec=((0:Nsc-1)*(1/Tc))';
c=PNSeqGen();
%% Generate user and MAI data
bits=round(rand(1,2*NoSymbs));
A= QPSKMod(bits,1,deg2rad(43)); 
ai1= Demux(A,width(A),Nsc);
MAI=zeros(M-1,length(A));
MAIres=cell(1,M-1);
rng default
for user=2:M
    bitsMAI= round(rand(1,2*NoSymbs));
    MAI(user-1,:)= QPSKMod(bitsMAI,0.2,deg2rad(10));
    MAIres{user-1}=Demux(MAI(user-1,:),width(A),Nsc);
end
ausers=[ai1, cell2mat(MAIres)];
ausers= unitymag(ausers); %every element with unity magnitude
ausers=ausers- real(mean(mean(ausers)));
%% 
[r,r_bar]=TxRxArr(lightvel,Fc,9);
%[r,r_bar]=TxRxArr(lightvel,Fc,"50");
N=length(r);
%[delays,beta,DODs,DOAs,VDops]= Channel_Param_Gen(1,0,Nsc);
[delays,beta,DODs,DOAs,VDops]=ChannelParam(0,0,0,Nsc);
% find f and gamma
f1j= computef(NoSymbs/Nsc,VDops(1,:),Fjvec,Fc,Tcs,lightvel,K);
f2j= computef(NoSymbs/Nsc,VDops(2,:),Fjvec,Fc,Tcs,lightvel,K);
f3j= computef(NoSymbs/Nsc,VDops(3,:),Fjvec,Fc,Tcs,lightvel,K);
f4j= computef(NoSymbs/Nsc,VDops(4,:),Fjvec,Fc,Tcs,lightvel,K);
f5j= computef(NoSymbs/Nsc,VDops(5,:),Fjvec,Fc,Tcs,lightvel,K);

f=[f1j,f2j,f3j,f4j,f5j];
%clear f1j f2j f3j f4j f5j;
gamma1= computegamma(beta(:,1:Nsc),DODs(1,:),Fjvec,r_bar,K);
gamma2= computegamma(beta(:,Nsc+1:2*Nsc),DODs(2,:),Fjvec,r_bar,K);
gamma3= computegamma(beta(:,2*Nsc+1:3*Nsc),DODs(3,:),Fjvec,r_bar,K);
gamma4= computegamma(beta(:,3*Nsc+1:4*Nsc),DODs(4,:),Fjvec,r_bar,K);
gamma5= computegamma(beta(:,4*Nsc+1:5*Nsc),DODs(5,:),Fjvec,r_bar,K);
gamma=[gamma1, gamma2,gamma3,gamma4,gamma5];
G= computeG(gamma,K);
%clear gamma1 gamma2 gamma3 gamma4 gamma5; 

SNR_abs= 10^(20/10);
PTx= 1;
Pnoise= abs( PTx/SNR_abs);
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
clear h1 h2 h3;
%% find x-eqn 17
x=zeros(N*Next,L);
for n=1:L
   x(:,n)= findX(ausers,f,gamma,H,J,M,Nsc,Nc,N,Next,K,n,L,Pnoise);
end
% load("x.mat");
noise= sqrt(Pnoise/2)* (randn(size(x))+1i*randn(size(x)));
x=x+noise;
%% eqn 24
[Rxx_theor,Rdes,Rmai,Risi,Rnn]= covtheor(H,G,N,Nc,K,Nsc,M,Pnoise); %compute theoretical cov. matx#
Rxx_prac= (1/L)* (x) * (x)';
%% DOA cost function
[Pn,~]= findPn(Rxx_theor,M*Nsc);
del_est=delays(1,:);
vel_est=VDops(1,:);
%[cost1d]=OneDCost(del_est,uk_est,Fjvec,r,Pn,J,c(:,1),Nsc,K);
[cost1d,DOAest]=experiment(del_est,vel_est, (1:360), Fjvec,r,Pn,J,c(:,1),Nsc,K);
Colours = {'red','g','blue'};
figure;
for k=1:K
    txt = ['Multipath ',num2str(k) ];
    plot(20*log10(cost1d(k,:)),'Color',Colours{k},'DisplayName',txt);
    hold on;
end
xlabel('DOA(degrees)'); ylabel('Gain(dB)'); title('DOA Estimation'); legend('show');
hold off;
%% Doppler STAR subspace and RAKE weights
j=1;
Hj=H(1:2*N*Nc*Nsc, (j-1)*K+1: j*K);
Gj= G(:, (j-1)*K+1: j*K);
gammaj= gamma(:,j);
Rdes=Hj*Gj*ctranspose(Hj);
out=findWeights(Rxx_prac,Hj,Gj,gammaj,K,N,Nc,Nsc);
weights= unitymag(out);
doppstar_weights= sum(weights,2,'omitnan');
%% 
SNIRout_doppstar= (ctranspose(doppstar_weights)*Rdes* (doppstar_weights))/ (ctranspose(doppstar_weights)*(Rnn)* (doppstar_weights));
10*log10(SNIRout_doppstar)
