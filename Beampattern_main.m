addpath('.\Utilities\');
addpath('.\Compute\');
%Define constants
NoSymbs=200;
M=5; %# users
N_bar=16;
Nsc=5;
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
% load("bitstream.mat");
% bits=bitstream(1,:);
bits= round(rand(1,2*NoSymbs));
phi_rad=deg2rad(43);
A= QPSKMod(bits,sqrt(2),phi_rad);
MAI= zeros(M-1,width(A));
%bitsMAI= bitstream(2,:);
bitsMAI= round(rand(1,2*NoSymbs));
MAI(1,:)= QPSKMod(bitsMAI,0.2,deg2rad(0)); %symbol streams have the same power

for user=2:M-1
    MAI(user,:)= shuffle(MAI(1,:));
end
ai1=Demux(A,width(A),Nsc);
MAIres=cell(1,M-1);
for user=2:M
    MAIres{user-1}=Demux(MAI(user-1,:),width(A),Nsc);
end
ausers=[ai1, cell2mat(MAIres)];
ausers= unitymag(ausers); %every element with unity magnitude
%make the signal zero mean
ausers=ausers- real(mean(mean(ausers)));
%% Define channel parameters
[r,r_bar]=TxRxArr(lightvel,Fc,"default");
%[r,r_bar]=TxRxArr(lightvel,Fc,"50");
N=length(r);
[delays,beta,DODs,DOAs,VDops]= Channel_Param_Gen(1,0);
% find f and gamma
f1j= computef(NoSymbs/Nsc,VDops(1,:),Fjvec,Fc,Tcs,lightvel,K);
f2j= computef(NoSymbs/Nsc,VDops(2,:),Fjvec,Fc,Tcs,lightvel,K);
f3j= computef(NoSymbs/Nsc,VDops(3,:),Fjvec,Fc,Tcs,lightvel,K);
f4j= computef(NoSymbs/Nsc,VDops(4,:),Fjvec,Fc,Tcs,lightvel,K);
f5j= computef(NoSymbs/Nsc,VDops(5,:),Fjvec,Fc,Tcs,lightvel,K);

f=[f1j,f2j,f3j,f4j,f5j];
clear f1j f2j f3j f4j f5j;
gamma1= computegamma(beta(:,1:5),DODs(1,:),Fjvec,r_bar,K);
gamma2= computegamma(beta(:,6:10),DODs(2,:),Fjvec,r_bar,K);
gamma3= computegamma(beta(:,11:15),DODs(3,:),Fjvec,r_bar,K);
gamma4= computegamma(beta(:,16:20),DODs(4,:),Fjvec,r_bar,K);
gamma5= computegamma(beta(:,21:25),DODs(5,:),Fjvec,r_bar,K);
gamma=[gamma1, gamma2,gamma3,gamma4,gamma5];
clear gamma1 gamma2 gamma3 gamma4 gamma5; 

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
%% find x-eqn 17
x=zeros(N*Next,L);
for n=1:L
   x(:,n)= findX(ausers,f,gamma,H,J,M,Nsc,Nc,N,Next,K,n,L,Pnoise);
end
% load("x.mat");
noise= sqrt(Pnoise/2)* (randn(size(x))+1i*randn(size(x)));
x=x+noise;
%% eqn 24
G= computeG(gamma,K);
Rxx_theor= covtheor(H,G,J,N,Nc,Nsc,M,Pnoise); %compute theoretical cov. matx
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
%[cost2d,del_est,uk_est]=faster2dcost(K,Nc,Nsc,delays(1,:),(0:140),akj,Fkj,Pn_res,J);
figure;
rv_range=(1:140);
delay_range=(0:Nc*Nsc-1);
surf(rv_range,delay_range,20*log10((cost2d)),'FaceAlpha',1,'EdgeAlpha',0.5);
xlabel('Velocity(m/s)'); ylabel('Delay(Ts s)'); zlabel('Gain(dB)');
title('Joint Delay-Doppler Velocity Estimation');
toc;
%% DOA cost function
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
%% Spatiotemporal beamformer weights
tic;
w=zeros(2*N*Nc*Nsc,Nsc);
for j=1:Nsc
    Hj= H(1:2*N*Nc*Nsc, (j-1)*K+1: j*K);
    Gj= G(:, (j-1)*K+1: j*K);
    gammaj= gamma(:,j);
    w(:,j)= subspaceWeightsj(Rxx_prac,Hj,Gj,gammaj,M,Nsc,Nc,N);
end
toc;
wjtotal= sum(w,2);
%% Doppler STAR vectors 
tic;
hdoppstarvecs= zeros(Nc*Nsc*2*N*Nc*Nsc,360*Nsc);
for delk= 0:Nc*Nsc-1
    for DOA=1:360
        for j=1:Nsc
            hdoppstar= DoppSTARmanifold(DOA,delk,VDops(1,1),J,Fjvec(j),Nsc,r,c(:,1));
            hdoppstarvecs( delk*2790+1: (delk+1)*2790, DOA+(j-1)*360)= hdoppstar; 
        end
    end
end
toc;
%% approach 1 ->use 1 subcarrier
gain=zeros(Nc*Nsc,360);
for delk=0:Nc*Nsc-1
    for DOA=1:360
        gain(delk+1,DOA)= ctranspose(w(:,1))* hdoppstarvecs( delk*2790+1: (delk+1)*2790, DOA);
    end
end
figure;
gain= abs(gain); 
surf((1:360),(0:Nc*Nsc-1),gain,'FaceAlpha',1,'EdgeAlpha',0.5);
shading('interp');
colormap('jet');
%% approach 1 -> sum across subcarriers
%h_overall= zeros(2790*Nc*Nsc,360);
for delk=0:Nc*Nsc-1 
    for DOA=1:360
        delcurr= hdoppstarvecs(delk*2790+1: (delk+1)*2790, DOA:360:end); %for a given delay,theta and Fj
        deldoa_sum= sum(delcurr,2); %sum across subcarriers to obtain a 2790x1 vector for delay,theta
       % h_overall(delk*2790+1: (delk+1)*2790,DOA)= deldoa_sum;
        gain(delk+1,DOA)= ctranspose(wjtotal)* deldoa_sum;
       % gain(delk+1,DOA)= ctranspose(w_RAKE_total)* deldoa_sum;
    end
end

figure;
gain= abs(gain); 
surf((1:360),(0:Nc*Nsc-1),gain,'FaceAlpha',1,'EdgeAlpha',0.5);
shading('interp');
colormap('jet');
%% space-only manifold beampattern ->figure 7b
%weights are 1s
%plotted with manifold as h=kron(S, J^delay * zero-padded user 1 pn code)
space_gain=[];
for delk=0:Nc*Nsc-1
    for DOA=1:360
        temp= getSpatiotemporalManifold(delk,DOA,J,c(:,1),r,Fc,0,lightvel);
        space_gain(delk+1,DOA)= wjtotal' * temp;
    end
end
figure;
surf(abs(space_gain));
shading('interp');
colormap('jet');

%alternative manifold -> S
%weights 1s
% for delk=0:Nc*Nsc-1
%     % for DOA=1:360
% %     manifold= computeManifoldRx(DOA,0,r,Fc,Fjvec(1),lightvel,"hwl");
% %     end
% space_gain(delk+1,DOA)= ones(1,2*N*Nc*Nsc) * manifold;
% end

%% Estimate of transmitted symbol vector
a_estim= w'*x(:,1);
%% obtaining equations 47-49
[w_RAKE,w_dec]=weights4749(H,gamma,K,N,Nc,Nsc);
w_RAKE_total= sum(w_RAKE,2);

%% new idea space
DOAest=[60,280,280];
space_gain= space_only_beampattern(DOAest,r,Fc,Fjvec,lightvel,Nc,Nsc);





