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
c=PNSeqGen();
bits= round(rand(1,2*NoSymbs)); %bitstream(1,:); % ;
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
[r,r_bar]=TxRxArr(lightvel,Fc);
[delays,beta,DODs,DOAs,VDops]= Channel_Param_Gen(0,0);
%find f and gamma
f1j= computef(ai1,VDops(1,:),Fjvec,Fc,Tcs,lightvel,K);
f2j= computef(ai1,VDops(2,:),Fjvec,Fc,Tcs,lightvel,K);
f3j= computef(ai1,VDops(3,:),Fjvec,Fc,Tcs,lightvel,K);
f4j= computef(ai1,VDops(4,:),Fjvec,Fc,Tcs,lightvel,K);

f=[f1j,f2j,f3j,f4j];
gamma1= computegamma(beta(:,1:5),DODs(1,:),Fjvec,r_bar,K);
gamma2= computegamma(beta(:,6:10),DODs(2,:),Fjvec,r_bar,K);
gamma3= computegamma(beta(:,11:15),DODs(3,:),Fjvec,r_bar,K);
gamma4= computegamma(beta(:,16:20),DODs(4,:),Fjvec,r_bar,K);
gamma=[gamma1, gamma2,gamma3,gamma4];
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
%% define MC simulation parameters
SNR=[0.5,3.1623,5,40,50,250,500,2500,5000];
L=200;
xaxis= L*SNR;
SNR_db=10.*log10(SNR);
numtrials=100;
RMSEgamma= zeros(numtrials,length(SNR));
RMSEDOA= zeros(numtrials,length(SNR));
%load("x.mat");
%% obtain noiseless x
x=zeros(N*Next,L);
for n=1:L
       x(:,n)=findX(ausers,f,gamma,H,J,M,Nsc,Nc,N,Next,K,n,L,0.1);
end
%% obtain x for every trial and SNR 
x_noise= cell(numtrials,length(SNR));
Pnoise= 1./SNR;
for trial=1:numtrials
    for ind=1:length(SNR)
        rng shuffle
        noise= sqrt(Pnoise(ind)/2)* (randn(size(x))+1i*randn(size(x)));
        x_noise{trial,ind}=x+noise;
    end
end
 %% obtain uk estimates for every SNR
 tic;
 RMSEradvel=zeros(numtrials,length(SNR));
%  rad_vel_range=(1:140);

 vel_est=[];
 
 for trial=1:numtrials
     uk_est=zeros(length(SNR),K);
     for ind=1:length(SNR)
        x_res= reshape(x_noise{trial,ind},2*Nc*Nsc,[]);
        Rxx_res= (1/width(x_res))* (x_res)*ctranspose(x_res);
        [Pn_res,~]= findPn(Rxx_res,length(Rxx_res)-M*Nsc);
        rad_vel_range=(1:140);
        [akj,Fkj]=findvecs(Fjvec,rad_vel_range,c(:,1),Nc,Nsc,Ts);
        [cost2d,del_est,uk_est(ind,:)]= faster2dcost(K,Nc,Nsc,delays(1,:),rad_vel_range,akj,Fkj,Pn_res,J);
        currRMSE=findRMSE(VDops(1,:),uk_est(ind,:));
        
        first=1;
        while(currRMSE==0)
            rad_vel_range=[];
            if(first==1)
                init=0.1;
                first=2;
            else
                init=init/10;
            end
            for k=1:K
                rad_vel_range=[rad_vel_range,uk_est(ind,k)-2:init:uk_est(ind,k)+2];
            end

            [akj,Fkj]=findvecs(Fjvec,rad_vel_range,c(:,1),Nc,Nsc,Ts);
            [cost2d,del_est,uk_est(ind,:)]= faster2dcost(K,Nc,Nsc,delays(1,:),rad_vel_range,akj,Fkj,Pn_res,J);
            currRMSE=findRMSE(VDops(1,:),uk_est(ind,:));
        end
        RMSEradvel(trial,ind)=currRMSE;
        
     end
     vel_est=[vel_est,uk_est];
 end
 toc;
%  %% precompute Pn matrices
%  Pn_mat=cell(numtrials,length(SNR));
%  for trial=1:numtrials
%      for ind=1:length(SNR)
%         xcurr= x_noise{trial,ind};
%         Rxx_prac= (1/L)* (xcurr) * ctranspose(xcurr);
%         [Pn,lambda_min]= findPn(Rxx_prac,M*Nsc);
%         Pn_mat{trial,ind}= Pn;
%      end
% 
% 
%  end
 %% obtain DOA estimates for every SNR
tic;
%load("vels100trials.mat");
init=0.01;
angles=[60,200,280]; %predicted values
theta=[];
for k=1:K
       theta=[theta, DOAs(1,k)-5:init: DOAs(1,k)+5];
end

 for trial=1:1
    for ind=1:3
        del_est=[140,110,30];
        %uk_est= [20,60,120];
        uk_est= vel_est(ind, (trial-1)*K+1: trial*K);
        xcurr= x_noise{trial,ind};
%         %use Rxx_theor
%         Rxx_theor= covtheor(H,G,J,N,Nc,Nsc,M,Pnoise(ind));
%[Pn,~]= findPn(Rxx_theor,M*Nsc);
        Rxx_prac= (1/L)* (xcurr) * ctranspose(xcurr);
        [Pn,lambda_min]= findPn(Rxx_prac,M*Nsc);
       

        %could use fmin search at this point
        [cost1d,~]=experiment(del_est,uk_est,theta,Fjvec,r,Pn,J,c(:,1),Nsc,K);
        DOAest= findMaxofPath(cost1d);
        DOAest= theta(DOAest);
        currerr= findRMSE(DOAs(1,:),DOAest);
        first=1;
        while(currerr==0)
            disp('zero');
            if(first==1)
                stepsize=init/10;
                first=2;
            else
                stepsize= stepsize/10;
            end
            %redefine search range
            thetas=[];
            for k=1:K
                thetas= [thetas, DOAest(k)-1: stepsize: DOAest(k)+1];
            end

            [cost1d,~]=experiment(del_est,uk_est,theta,Fjvec,r,Pn,J,c(:,1),Nsc,K);
            DOAest= findMaxofPath(cost1d);
            DOAest= theta(DOAest);
            currerr= findRMSE(DOAs(1,:),DOAest);
        end
        RMSEDOA(trial,ind)= currerr;
        
    end
    disp(trial);
 end
toc;

%% CRB 
% xaxis= SNR*L;
% dsmag= manifoldDer(60,0,r,Fc,Fjvec(1),lightvel,"hwl");
% CRB= 1./(2.*xaxis.*(dsmag^2));
% CRB_Dop= (1/Nc*Nsc)*CRB;
% figure;
% loglog(xaxis,CRB_Dop,'o','DisplayName','CRB');
% xlim([10e2,10e7]);
% % plot(CRB);
%% 
%the xaxis2 vector goes from 10^2 to 10^5
figure;
loglog(xaxis,mean(RMSEradvel(1:2,:)),'o','DisplayName','Radial Vel RMSE');
xlabel('SNRxL'); ylabel('Estimation RMSE');
grid on;
ylim([10e-5 10e1]); legend('show');
%    [cost2d,del_est,uk_est]= TwoDcost(K,Nc,Nsc,delays(1,:),akj,Fkj,Pn_res,J);
% %     uk_est(3)=118;
%     RMSEradvel(1,ind)= findRMSE(VDops(1,:),uk_est);

