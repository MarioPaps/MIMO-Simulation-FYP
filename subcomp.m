addpath('.\Utilities\');
addpath('.\Compute\');
addpath('.\Beampattern_files\');
M=3; %# users   
N_bar=16;
L=40; %number of snaphots of x[n]
Nc=31;
Ts=0.1e-6;
Fc=20e9;
K=3; %3 paths per user
lightvel=3e8;
lambda= lightvel/Fc;
c=PNSeqGen();
%% trial code
numtrials=100;
% SNIR_nsc=zeros(numtrials,10);
% up=[];
% down=[];
for trial=1:1
    for Nsc=9:10
        Tc=Ts*Nsc;
        Tcs=Nc*Tc; 
        Fjvec=((0:Nsc-1)*(1/Tc))';
        NoSymbs= L*Nsc;
        %generate signal
        bits=round(rand(1,2*NoSymbs));
        A= QPSKMod(bits,1,deg2rad(43)); 
        ai1= Demux(A,width(A),Nsc);
        MAI=zeros(M-1,length(A));
        MAIres=cell(1,M-1);
        rng default
        for user=2:M
            bitsMAI= round(rand(1,2*NoSymbs));
            MAI(user-1,:)= QPSKMod(bitsMAI,0.002,deg2rad(10));
            MAIres{user-1}=Demux(MAI(user-1,:),width(A),Nsc);
        end
        ausers=[ai1, cell2mat(MAIres)];
        ausers= unitymag(ausers); %every element with unity magnitude
        ausers=ausers- real(mean(mean(ausers)));
    
        [r,r_bar]=TxRxArr(lightvel,Fc,20);
        N=length(r);
    
        [delays,beta,DODs,DOAs,VDops]=ChannelParam(0,0,1,Nsc);
        
        f1j= computef(NoSymbs/Nsc,VDops(1,:),Fjvec,Fc,Tcs,lightvel,K);
        f2j= computef(NoSymbs/Nsc,VDops(2,:),Fjvec,Fc,Tcs,lightvel,K);
        f3j= computef(NoSymbs/Nsc,VDops(3,:),Fjvec,Fc,Tcs,lightvel,K);
      %  f4j= computef(NoSymbs/Nsc,VDops(4,:),Fjvec,Fc,Tcs,lightvel,K);
       % f5j= computef(NoSymbs/Nsc,VDops(5,:),Fjvec,Fc,Tcs,lightvel,K);
        f=[f1j,f2j,f3j];
        clear f1j f2j f3j ;
        gamma1= computegamma(beta(:,1:Nsc),DODs(1,:),Fjvec,r_bar,K);
        gamma2= computegamma(beta(:,Nsc+1:2*Nsc),DODs(2,:),Fjvec,r_bar,K);
        gamma3= computegamma(beta(:,2*Nsc+1:3*Nsc),DODs(3,:),Fjvec,r_bar,K);
     %   gamma4= computegamma(beta(:,3*Nsc+1:4*Nsc),DODs(4,:),Fjvec,r_bar,K);
    %    gamma5= computegamma(beta(:,4*Nsc+1:5*Nsc),DODs(5,:),Fjvec,r_bar,K);
        gamma=[gamma1, gamma2,gamma3];
        G= computeG(gamma,K);
        clear gamma1 gamma2 gamma3 ; 
    
        SNR_abs= 10^(20/10);
        PTx= 1;
        Pnoise= abs( PTx/SNR_abs);
    
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
    
        x=zeros(N*Next,L);
        for n=1:L
           x(:,n)= findX(ausers,f,gamma,H,J,M,Nsc,Nc,N,Next,K,n,L,Pnoise);
        end
        % load("x.mat");
        noise= sqrt(Pnoise/2)* (randn(size(x))+1i*randn(size(x)));
        x=x+noise;
        [Rxx_theor,Rdes,Rmai,Risi,Rnn]= covtheor(H,G,N,Nc,K,Nsc,M,Pnoise); %compute theoretical cov. matx#
        Rxx_prac= (1/L)* (x) * (x)';
        %Doppler weights
        j=1;
        Hj=H(1:2*N*Nc*Nsc, (j-1)*K+1: j*K);
        Gj= G(:, (j-1)*K+1: j*K);
        gammaj= gamma(:,j);
        Rdes=Hj*Gj*ctranspose(Hj);
        out=findWeights(Rxx_prac,Hj,Gj,gammaj,K,N,Nc,Nsc);
        weights= unitymag(out);
        doppstar_weights= sum(weights,2,'omitnan');
        up(trial,Nsc)=ctranspose(doppstar_weights)*Rdes* (doppstar_weights);
        down(trial,Nsc)=ctranspose(doppstar_weights)*(Rxx_prac-Rdes)* (doppstar_weights);

        SNIRout_doppstar= (ctranspose(doppstar_weights)*Rdes* (doppstar_weights))/ (ctranspose(doppstar_weights)*(Rnn)* (doppstar_weights));
        SNIR_nsc(trial,Nsc)=10*log10(real(SNIRout_doppstar));
       
    end
end
%% 
 plot((1:10),10*log10(real(up(1,:)./down(1,:))),'-o');
 figure;
 plot((1:10),10*log10(real(up(1,:))),'-o');
 figure;
 plot((1:10),10*log10(real(down(1,:))),'-o');
 %% 
 SNIR_sub= up;
 SNIR_sub(2,9)=SNIR_sub(1,9);
 SNIR_sub(2,10)=SNIR_sub(1,10);
 save('SNIR_sub.mat','SNIR_sub');

