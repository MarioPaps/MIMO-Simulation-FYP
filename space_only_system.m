addpath('.\Utilities\');
addpath('.\Compute\');
addpath('.\Beampattern_files\');
%Define constants
NoSymbs=200;
M=2; %# users
N_bar=16;
Nsc=1;
%L=NoSymbs/Nsc;
L=40;
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
%%  Generate user and MAI data
NFR= 10.^((0:10:60)./10); %NFR levels
bits=round(rand(1,2*NoSymbs));
A= QPSKMod(bits,1,deg2rad(43)); 
MAI= MAIfromNFR(A,NFR(7),M,NoSymbs);
ai1= Demux(A,width(A),Nsc);
MAIres=cell(1,M-1);
for user=2:M
    MAIres{user-1}=Demux(MAI(user-1,:),width(A),Nsc);
end
ausers=[ai1, cell2mat(MAIres)];
ausers= unitymag(ausers); %every element with unity magnitude
ausers=ausers- real(mean(mean(ausers)));
%% channel parameters
[r,r_bar]=TxRxArr(lightvel,Fc,100);
N=length(r);
[delays,beta,DODs,DOAs,VDops]= Channel_Param_Gen(0,0,Nsc);
% find gamma
gamma1= computegamma(beta(:,1:5),DODs(1,:),Fjvec,r_bar,K);
gamma2= computegamma(beta(:,6:10),DODs(2,:),Fjvec,r_bar,K);
gamma3= computegamma(beta(:,11:15),DODs(3,:),Fjvec,r_bar,K);
gamma4= computegamma(beta(:,16:20),DODs(4,:),Fjvec,r_bar,K);
gamma5= computegamma(beta(:,21:25),DODs(5,:),Fjvec,r_bar,K);
gamma=[gamma1, gamma2,gamma3,gamma4,gamma5];
% clear gamma1 gamma2 gamma3 gamma4 gamma5; 

SNR_abs= 10^(20/10);
Pnoise= abs(1/SNR_abs);
%% space-only received signal
x=zeros(N,L);
for n=1:L
    x(:,n)=findXspace(ausers,gamma,DOAs,r,Fjvec,n,L,M,K,Fc);
end
x=x+sqrt(Pnoise/2)* (randn(size(x))+1i*randn(size(x)));
Rxx= (1/width(x))*x*ctranspose(x);
%% plot MUSIC spectrum to check correctness
musicout= MUSIC_cost(Rxx,r,(1:360),Fjvec,M,Fc);
figure;
plot((1:360),20*log10(abs(musicout)));
%% weights computation
space_weights=findSubWeightsSpace(Rxx,DOAs(1,:),gamma,r,Fjvec,Fc);
space_weights=unitymag(space_weights);
space_weights=sum(space_weights,2);
%% Rjj computation
amai= ausers(:,width(A)+1:end);
amai= reshape(amai,[],M-1);
correlations= corrcoef(amai(1:end));

%% nfr computation
Sdd= computeManifoldRx(DOAs(1,:)',0,r,Fc,Fjvec(1),lightvel);
Rdd= Sdd * Sdd';
SNIRout= (space_weights'* Rdd * space_weights)/(space_weights'*(Rxx-Rdd)*space_weights)

