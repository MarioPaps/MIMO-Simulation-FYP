%find SNIRout in terms of total number of paths for all users
addpath('.\Utilities\');
addpath('.\Compute\');
%Define constants
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
%% transmitted signal generation
bits= round(rand(1,2*NoSymbs));
A= QPSKMod(bits,sqrt(2),deg2rad(43));
ai1=Demux(A,width(A),Nsc);

MAI=zeros(M-1,length(A));
MAIres=cell(1,M-1);
rng default
for user=2:M
    bitsMAI= round(rand(1,2*NoSymbs));
    MAI(user-1,:)= QPSKMod(bitsMAI,1,deg2rad(10));
    MAIres{user-1}=Demux(MAI(user-1,:),width(A),Nsc);
end

ausers=[ai1, cell2mat(MAIres)];
ausers= unitymag(ausers); %every element with unity magnitude
ausers=ausers- real(mean(mean(ausers)));
%% Obtain SNIRout for varying K
[r,r_bar]=TxRxArr(lightvel,Fc,5);
N=length(r);
K=[3:3:60];
numtrials=10;
%results4= zeros(numtrials,length(K));
for trial=11:50
    count=1;
    for k=K
        betaprime=[0.5,0.1,0.1,0.1,0.1];
        beta= repmat(betaprime,k);
        delays=randi([2,Nc],M,k);
        DODs=randi([1,360],M,k);
        DOAs=randi([1,360],M,k);
        VDops=randi([1,140],M,k);
        f=zeros(k*Nsc,L*M);
        gamma=zeros(k,M*Nsc);
        for user=1:M
            f(:,(user-1)*L+1: user*L)=computef(NoSymbs/Nsc,VDops(user,:),Fjvec,Fc,Tcs,lightvel,k);
            gamma(:,(user-1)*Nsc+1:user*Nsc)=computegamma(beta(:,(user-1)*Nsc+1: user*Nsc),DODs(user,:),Fjvec,r_bar,k);
        end
        G= computeG(gamma,k); 
        SNR_abs= 10^(20/10);
        PTx= 1;
        Pnoise= abs( PTx/SNR_abs);
        J = [zeros(1,2*Nc*Nsc-1) 0; eye(2*Nc*Nsc-1), zeros(2*Nc*Nsc-1,1)];
        Next= 2*Nc*Nsc;
        H=zeros(M*2*N*Nc*Nsc,k*Nsc);
        for i=1:M
            for j=1:Nsc
                row_start= (i-1)*N*Next+1;
                row_end= i*N*Next;
                for kter=1:k
                      h= DoppSTARmanifold(DOAs(i,kter),delays(i,kter),VDops(i,kter),J,Fjvec(j),Nsc,r,c(:,i));
                      H(row_start:row_end,(j-1)*k+kter)=h;
                end
            end
        end
    %received signal
    x=zeros(N*Next,L);
    for n=1:L
        x(:,n)=findX(ausers,f,gamma,H,J,M,Nsc,Nc,N,Next,k,n,L,Pnoise);
    end
    rng shuffle;
    noise= sqrt(Pnoise/2)* (randn(size(x))+1i*randn(size(x)));
    x=x+noise;
    Rxx_prac= (1/L)* (x) * (x)';
    %[~,~,Rmai,Risi,~]=covtheor(H,G,N,Nc,k,Nsc,M,Pnoise);
    Rnn=Pnoise*speye(2*N*Nc*Nsc);
    j=1;
    Hj=H(1:2*N*Nc*Nsc, (j-1)*k+1: j*k);
    Gj= G(:, (j-1)*k+1: j*k);
    gammaj= gamma(:,j);
    Rdes=Hj*Gj*ctranspose(Hj);
    out=findWeights(Rxx_prac,Hj,Gj,gammaj,k,N,Nc,Nsc);
    weights= unitymag(out);
    doppstar_weights= sum(weights,2,'omitnan');
    SNIRout_doppstar= (ctranspose(doppstar_weights)*Rdes* (doppstar_weights))/ (ctranspose(doppstar_weights)*(Rnn)* (doppstar_weights));
    results4(trial,count)=10*log10(SNIRout_doppstar);
    count=count+1;
    end

end
%% 
overall=[results;results2;results3];
%overall=results3;
figure;
plot(K,mean(abs(overall)),'-o');
grid on;
grid minor;
ylim([-50 50]);
%% 
figure;
plot(K,aa,'-o');
grid on;
grid minor;
ylim([-50 50]);
title('Effect of Total Number of Paths on SNIRout');
xlabel('Paths Per User'); ylabel('SNIRout (dB)');
%% 
figure;
plot(K,mean(real(results4)),'-o');
grid on;
grid minor;
ylim([-50 50]);
title('Effect of Total Number of Paths on SNIRout');
xlabel('Paths Per User'); ylabel('SNIRout (dB)');
