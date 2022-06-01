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
%% Define channel parameters
[r,r_bar]=TxRxArr(lightvel,Fc,9);
N=length(r);
[delays,beta,DODs,DOAs,VDops]= ChannelParam(0,0,0,Nsc);
% find f and gamma
f1j= computef(L,VDops(1,:),Fjvec,Fc,Tcs,lightvel,K);
f2j= computef(L,VDops(2,:),Fjvec,Fc,Tcs,lightvel,K);
f3j= computef(L,VDops(3,:),Fjvec,Fc,Tcs,lightvel,K);
f4j= computef(L,VDops(4,:),Fjvec,Fc,Tcs,lightvel,K);
f5j= computef(L,VDops(5,:),Fjvec,Fc,Tcs,lightvel,K);
f=[f1j,f2j,f3j,f4j,f5j];
clear f1j f2j f3j f4j f5j;

gamma1= computegamma(beta(:,1:Nsc),DODs(1,:),Fjvec,r_bar,K);
gamma2= computegamma(beta(:,Nsc+1:2*Nsc),DODs(2,:),Fjvec,r_bar,K);
gamma3= computegamma(beta(:,2*Nsc+1:3*Nsc),DODs(3,:),Fjvec,r_bar,K);
gamma4= computegamma(beta(:,3*Nsc+1:4*Nsc),DODs(4,:),Fjvec,r_bar,K);
gamma5= computegamma(beta(:,4*Nsc+1:5*Nsc),DODs(5,:),Fjvec,r_bar,K);
gamma=[gamma1, gamma2,gamma3,gamma4,gamma5];
G= computeG(gamma,K);
clear gamma1 gamma2 gamma3 gamma4 gamma5; 
SNR_abs= 10^(20/10);
PTx= 1;
Pnoise= abs( PTx/SNR_abs);
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
%         for k=1:K
%             h= DoppSTARmanifold(DOAs(i,k),delays(i,k),VDops(i,k),J,Fjvec(j),Nsc,r,c(:,i));
%         end
        h1= DoppSTARmanifold(DOAs(i,1),delays(i,1),VDops(i,1),J,Fjvec(j),Nsc,r,c(:,i));
        h2= DoppSTARmanifold(DOAs(i,2),delays(i,2),VDops(i,2),J,Fjvec(j),Nsc,r,c(:,i));
        h3= DoppSTARmanifold(DOAs(i,3),delays(i,3),VDops(i,3),J,Fjvec(j),Nsc,r,c(:,i));
        H(row_start:row_end,col_start:col_end)=[h1,h2,h3];
    end
end
clear h1 h2 h3;
%% user 1 transmitted signal
bits=round(rand(1,2*NoSymbs));
A= QPSKMod(bits,1,deg2rad(43)); 
ai1= Demux(A,width(A),Nsc);
%% MC simulation for Doppler STAR subspace beamformer
numtrials=200;
NFR= 10.^((0:10:60)./10); %NFR levels
doppstar_result=zeros(numtrials,length(NFR));
for trial=1:100
  [~,MAInfr,MAIpowers]=MAIfromNFR(A,NFR,M,NoSymbs);
    beta_nfr= betas_NFR(MAIpowers,M,Nsc,K);
    for ind=1:length(NFR)
        holder=zeros(Nsc,(NoSymbs/Nsc)*(M-1));
        for it=1:M-1
           holder(:,(it-1)*(NoSymbs/Nsc)+1: it*(NoSymbs/Nsc))=Demux(MAInfr{1,ind}(it,:),width(A),Nsc);
        end

        ausers=[ai1,holder];
        ausers= unitymag(ausers); 
       
        gamma1= computegamma(beta(:,1:Nsc),DODs(1,:),Fjvec,r_bar,K);
        gamma2= computegamma(beta(:,Nsc+1:2*Nsc),DODs(2,:),Fjvec,r_bar,K);
        gamma3= computegamma(beta(:,2*Nsc+1:3*Nsc),DODs(3,:),Fjvec,r_bar,K);
        gamma4= computegamma(beta(:,3*Nsc+1:4*Nsc),DODs(4,:),Fjvec,r_bar,K);
        gamma5= computegamma(beta(:,4*Nsc+1:5*Nsc),DODs(5,:),Fjvec,r_bar,K);
        gamma=[gamma1, gamma2,gamma3,gamma4,gamma5];
        clear gamma1 gamma2 gamma3 gamma4 gamma5;

        %received signal
        x=zeros(N*Next,L);
        for n=1:L
            x(:,n)=findX(ausers,f,gamma,H,J,M,Nsc,Nc,N,Next,K,n,L,Pnoise);
        end
        % load("x.mat");
        rng shuffle;
        noise= sqrt(Pnoise/2)* (randn(size(x))+1i*randn(size(x)));
        x=x+noise;
            
        %covariance matrix
        [Rxx_theor,Rdes,Rmai,Risi,Rnn]= covtheor(H,G,N,Nc,K,Nsc,M,Pnoise); %compute theoretical cov. matx#
        Rxx_prac= (1/L)* (x) * (x)';
    
        j=1;
        Hj=H(1:2*N*Nc*Nsc, (j-1)*K+1: j*K);
        Gj= G(:, (j-1)*K+1: j*K);
        gammaj= gamma(:,j);
        Rdes=Hj*Gj*ctranspose(Hj);
        out=findWeights(Rxx_prac,Hj,Gj,gammaj,K,N,Nc,Nsc);
        weights= unitymag(out);
        doppstar_weights= sum(weights,2,'omitnan');
        
        % SNIRout evaluation
        SNIRout_doppstar= (ctranspose(doppstar_weights)*Rdes* (doppstar_weights))/ (ctranspose(doppstar_weights)*(Rnn)* (doppstar_weights));
        doppstar_result(trial,ind)= 10*log10(real(SNIRout_doppstar));
        save('doppstar_result.mat','doppstar_result');
    end  
end

%% 
%semilogy(10*log10(NFR),mean(doppstar_result(1:100,:)),'-o');
%stem(mean(doppstar_result(1:100,:)));
%% MC simulation for Doppler STAR RAKE beamformer
numtrials=100;
NFR= 10.^((0:5:60)./10); %NFR levels
%rake_result=zeros(numtrials,length(NFR));

up=[];
downer=[];

for trial=4:100
    [~,MAInfr,MAIpowers]=MAIfromNFR(A,NFR,M,NoSymbs);
    beta_nfr= betas_NFR(MAIpowers,M,Nsc,K);
    
    for ind=1:length(NFR)
        holder=zeros(Nsc,(NoSymbs/Nsc)*(M-1));
        for it=1:M-1
           holder(:,(it-1)*(NoSymbs/Nsc)+1: it*(NoSymbs/Nsc))=Demux(MAInfr{1,ind}(it,:),width(A),Nsc);
        end
        ausers=[ai1,holder];
        ausers= unitymag(ausers); %every element with unity magnitude

        %modify channel parameters
        beta= cell2mat(beta_nfr(ind));

        gamma1= computegamma(beta(:,1:Nsc),DODs(1,:),Fjvec,r_bar,K);
        gamma2= computegamma(beta(:,Nsc+1:2*Nsc),DODs(2,:),Fjvec,r_bar,K);
        gamma3= computegamma(beta(:,2*Nsc+1:3*Nsc),DODs(3,:),Fjvec,r_bar,K);
        gamma4= computegamma(beta(:,3*Nsc+1:4*Nsc),DODs(4,:),Fjvec,r_bar,K);
        gamma5= computegamma(beta(:,4*Nsc+1:5*Nsc),DODs(5,:),Fjvec,r_bar,K);
        gamma=[gamma1, gamma2,gamma3,gamma4,gamma5];
        clear gamma1 gamma2 gamma3 gamma4 gamma5;
        G=computeG(gamma,K);

        %received signal
        x=zeros(N*2*Nc*Nsc,L);
        for n=1:L
            x(:,n)=findX(ausers,f,gamma,H,J,M,Nsc,Nc,N,Next,K,n,L,Pnoise);
        end
        rng shuffle;
        noise= sqrt(Pnoise/2)* (randn(size(x))+1i*randn(size(x)));
        x=x+noise;
        %DOA estimation    
        [Rxx_theor,Rdes,Rmai,Risi,Rnn]=covtheor(H,G,N,Nc,K,Nsc,M,Pnoise);
        [w_rake,~]= weights4749(H,gamma,K,N,Nc,Nsc);
%         up=[up,ctranspose(w_rake)*Rdes* (w_rake)];
%         downer=[downer,ctranspose(w_rake)*(Rmai+Risi+Rnn)* (w_rake)];

        SNIRout_rake= (ctranspose(w_rake)*Rdes* (w_rake))/ (ctranspose(w_rake)*(Rmai+Risi+Rnn)* (w_rake));
        rake_result(trial,ind)= 10*log10(real(SNIRout_rake));
        
    end
end
%% 
plot(10*log10(NFR),mean(doppstar_result(1:100,:)),'-o'); 
hold on;
plot(10*log10(NFR),mean(rake_result(1,:)),'-o'); %add 25 dB worst case scenario
hold on;
plot(10*log10(NFR),mean(spacesub_result(1:50,:)),'-o')
%% space-only system with subspace weights
numtrials=100;
NFR= 10.^((0:5:60)./10); %NFR levels
%spacesub_result=zeros(numtrials,length(NFR));
space2=zeros(numtrials,length(NFR));
up=[];
downer=[];

for trial=1:10
    [~,MAInfr,MAIpowers]=MAIfromNFR(A,NFR,M,NoSymbs);
    beta_nfr= betas_NFR(MAIpowers,M,Nsc,K);
    for ind=1:length(NFR)
        holder=zeros(Nsc,(NoSymbs/Nsc)*(M-1));
        for it=1:M-1
           holder(:,(it-1)*(NoSymbs/Nsc)+1: it*(NoSymbs/Nsc))=Demux(MAInfr{1,ind}(it,:),width(A),Nsc);
        end
        ausers=[ai1,holder];
        ausers= unitymag(ausers); %every element with unity magnitude
       
        %generate received signal
        x=zeros(N,L);
        for n=1:L
            x(:,n)=findXspace(ausers,gamma,DOAs,r,Fjvec,n,L,M,K,Fc);
        end
        %rng shuffle;
        Pnoise=0.01;
        x=x+sqrt(Pnoise/2)* (randn(size(x))+1i*randn(size(x)));
        Rxx= (1/width(x))*x*ctranspose(x);
        
        %weights computation
        space_weights=findSubWeightsSpace(Rxx,DOAs(1,:),gamma,r,Fjvec,Fc);
        space_weights=unitymag(space_weights);
        space_weights=sum(space_weights,2);

        Rjj=spaceRjj(DOAs,r,Fjvec,Fc,M);

        Sdd= computeManifoldRx(DOAs(1,:)',0,r,Fc,Fjvec(1),lightvel);
        %Rdd= Sdd *gamma(1,1)* Sdd';
        Rdd= Sdd*Sdd';

        up=[up,space_weights'* Rdd * space_weights];
        downer=[downer,space_weights'*(Rxx-Rdd)*space_weights];
       % SNIRout= (space_weights'* Rdd * space_weights)/(space_weights'*(Rxx-Rdd)*space_weights);
        SNIRout= (space_weights'* Rdd * space_weights)/(space_weights'*(Rjj+0.01*speye(N,N))*space_weights);

        %spacesub_result(trial,ind)=10*log10(real(SNIRout));
        space2(trial,ind)= 10*log10(real(SNIRout_rake));
    end

end
%% 
semilogy(10*log10(NFR),mean(spacesub_result(1:10,:)),'-o')


