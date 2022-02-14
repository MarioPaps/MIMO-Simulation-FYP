addpath('.\Utilities\');
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
light_vel=3e8;

Fjvec=((0:Nsc-1)*(1/Tc))';
%% Generate user and MAI data
%create 400 bits long bitstream -> 200 channel symbols via QPSK modulation
load("bitstream.mat");
radius=sqrt(2);
phi_rad= 43*(pi/180);
numsymbs=4;
symbols= radius*exp(1i*(phi_rad+(0:numsymbs-1)*(pi/2)));
A= QPSKMod(bitstream(1,:),radius,phi_rad); %desired user's symbols
MAI= zeros(M-1,NoSymbs);
MAIsyms=[0.01+0.1i;0.01-0.1i];
%Generate MAI
for user=1:M-1
    id1=randi([1,100],1,100);
    id2=randi([101,200],1,100);
    
    MAI(user,id1)= MAIsyms(1);
    MAI(user,id2)= MAIsyms(2);
end
%demultiplexer output for every user
%The 200 symbols of every user are reshaped into a 5*col matrix
a_i1= Demux(A(1,:),M,NoSymbs,Nsc);
a_i2= Demux(MAI(1,:),M,NoSymbs,Nsc);
a_i3= Demux(MAI(2,:),M,NoSymbs,Nsc);
a_i4= Demux(MAI(3,:),M,NoSymbs,Nsc);
a_i5= Demux(MAI(4,:),M,NoSymbs,Nsc);
disp('demux out');
%gold codes
load("pnseq.mat");
%% antenna array setup
rx_bar=[2.56,2.37,1.81,0.98,0,-0.98,-1.81,-2.37,-2.56,-2.37,-1.81,-0.98,0,0.98,1.81,2.37]';
ry_bar=[0,0.98,1.81,2.37,2.56,2.37,1.81,0.98,0,-0.98,-1.81,-2.37,-2.56,-2.37,-1.81,-0.98]';
rz_bar=zeros(16,1);
r_bar=[rx_bar ry_bar rz_bar];
rx=[1.43 0.90 -0.05 -0.98 -1.45 -1.24 -0.45 0.55 1.29]';
ry=[0.30 1.15 1.46 1.09 0.20 -0.77 -1.39 -1.36 -0.69]';
rz=zeros(N,1);
r=[rx ry rz];
%% Define channel parameters
[delays,beta,DODs,DOAs,VDops]= Channel_Param_Gen();
beta= cell2mat(beta);
f_users=cell(1,M);
f1j= computef(a_i1,VDops(1,:),Fjvec,Fc,Tcs,light_vel,N_bar);
f2j= computef(a_i2,VDops(2,:),Fjvec,Fc,Tcs,light_vel,N_bar);
f3j= computef(a_i3,VDops(3,:),Fjvec,Fc,Tcs,light_vel,N_bar);
f4j= computef(a_i4,VDops(4,:),Fjvec,Fc,Tcs,light_vel,N_bar);
f5j= computef(a_i5,VDops(5,:),Fjvec,Fc,Tcs,light_vel,N_bar);

f=[f1j,f2j,f3j,f4j,f5j];
gamma1= computegamma(beta(:,1:5),DODs(1,:),Fjvec,r_bar);
gamma2= computegamma(beta(:,6:10),DODs(2,:),Fjvec,r_bar);
gamma3= computegamma(beta(:,11:15),DODs(3,:),Fjvec,r_bar);
gamma4= computegamma(beta(:,16:20),DODs(4,:),Fjvec,r_bar);
gamma5= computegamma(beta(:,21:25),DODs(5,:),Fjvec,r_bar);
gamma=[gamma1, gamma2,gamma3,gamma4,gamma5];

ausers=[a_i1,a_i2,a_i3,a_i4,a_i5];
%% equation 17
J = [zeros(1,2*Nc*Nsc-1) 0; eye(2*Nc*Nsc-1), zeros(2*Nc*Nsc-1,1)];
n=2;
upper=40;
Next= 2*Nc*Nsc;
x=zeros(2*N*Nc*Nsc,upper);
H=zeros(M*2*N*Nc*Nsc,K*Nsc);
for i=1:M
    for j=1:Nsc
        row_start= (i-1)*N*Next+1;
        row_end= i*N*Next;
        col_start= (j-1)*K+1;
        col_end= j*K;
        h1= DoppSTARmanifold(DOAs(i,1),delays(i,1),VDops(i,1),J,Fjvec(j),r,c(:,i));
        h2= DoppSTARmanifold(DOAs(i,2),delays(i,2),VDops(i,2),J,Fjvec(j),r,c(:,i));
        h3= DoppSTARmanifold(DOAs(i,3),delays(i,3),VDops(i,3),J,Fjvec(j),r,c(:,i));
        H(row_start:row_end,col_start:col_end)=[h1,h2,h3];
%         for path=1:K
%             hijk=DoppSTARmanifold(DOAs(i,path),delays(i,path),VDops(i,path),J,Fjvec(j),r,c(:,i));
%             H(row_start:row_end,col_start:col_end)=hijk;
%         end
    end
end

%% eqn17

for n=1:upper
 count=1;
 product=zeros(2*N*Nc*Nsc,M*Nsc);
for i=1:5
    for j=1:5
        row_start= (i-1)*N*Next+1;
        row_end= i*N*Next;
        col_start= (j-1)*K+1;
        col_end= j*K;
        %1st matrix
%         Hij=zeros(2*N*Nc*Nsc,K);
        Hij= H(row_start:row_end,col_start:col_end);
%         for path=1:K
%             Hij(:,path)=DoppSTARmanifold(DOAs(i,path),delays(i,path),VDops(i,path),J,Fjvec(j),r,c(:,i));
%         end
        firstmat=[Hij, computeHijcomb(Hij,J)];

        %2nd matrix
        if(n==1)
            currsymb= ausers(j,(i-1)*upper+n);
            nextsymb= ausers(j,(i-1)*upper+n+1);
            secondmat=[(gamma(:,(i-1)*M+j).*f((j-1)*K+1: j*K,(i-1)*upper+n))* currsymb;
                     zeros(3,1);
                   (gamma(:,(i-1)*M+j).*f((j-1)*K+1: j*K,(i-1)*upper+n+1))* nextsymb];
        elseif(n==upper)
            currsymb= ausers(j,(i-1)*upper+n);
            secondmat=[(gamma(:,(i-1)*M+j).*f((j-1)*K+1: j*K,(i-1)*upper+n))* currsymb;
                     (gamma(:,(i-1)*M+j).*f((j-1)*K+1: j*K,(i-1)*upper+n-1))* prevsymb;
                     zeros(3,1)];
        else
            currsymb= ausers(j,(i-1)*upper+n);
            prevsymb= ausers(j,(i-1)*upper+n-1);
            nextsymb= ausers(j,(i-1)*upper+n+1);
            secondmat=[(gamma(:,(i-1)*M+j).*f((j-1)*K+1: j*K,(i-1)*upper+n))* currsymb;
                       (gamma(:,(i-1)*M+j).*f((j-1)*K+1: j*K,(i-1)*upper+n-1))* prevsymb;
                       (gamma(:,(i-1)*M+j).*f((j-1)*K+1: j*K,(i-1)*upper+n+1))* nextsymb];
        end
        %multiplication
        product(:,count)= firstmat * secondmat;
        count=count+1;

%        

        
% %         second= gamma(:,(i-1)*M+j)
% %         third= f((j-1)*K+1: j*K,(i-1)*upper+n)
%         currsymb= ausers(j,(i-1)*upper+n);
% %         prevsymb= ausers(j,(i-1)*upper+n-1);
% %         nextsymb= ausers(j,(i-1)*upper+n+1);
%         secondmat=[(gamma(:,(i-1)*M+j).*f((j-1)*K+1: j*K,(i-1)*upper+n))* currsymb;
%                    (gamma(:,(i-1)*M+j).*f((j-1)*K+1: j*K,(i-1)*upper+n-1))* prevsymb;
%                    (gamma(:,(i-1)*M+j).*f((j-1)*K+1: j*K,(i-1)*upper+n+1))* nextsymb];
        %multiplication
%         product(:,count)= firstmat * secondmat;
%         count=count+1;
    end
end
x(:,n)= sum(product,2);
disp(n);
end
%% 
SNR_abs= 10^(20/10);
P_Tx= (1/length(bitstream(1,:)))*bitstream(1,:)* bitstream(1,:)';
Pnoise= abs( P_Tx/SNR_abs);
noise= sqrt(Pnoise/2)* (randn(size(x))+1i*randn(size(x)));
x_noisy= x+noise;
%% eqn 24
G= computeG(gamma,K);
Rxx_theor= covtheor(H,G,J,N,Nc,Nsc,Pnoise);
Rxx_prac= (1/width(x_noisy))* (x_noisy) * (x_noisy)';
%% 2d cost function
[akj,Fkj]=findvecs(Fjvec,c(:,1),Nc,Nsc,Ts);
x_res= reshape(x_noisy,2*Nc*Nsc,[]);
Rxx_res=(1/width(x_res))* (x_res)*ctranspose(x_res);
Pn_res= findPn(Rxx_res,M);
[cost2d,del_est,uk_est]= TwoDcost(N,Nc,Nsc,Fjvec,akj,Fkj,Pn_res,J);

figure;
surf(20*log10(abs(cost2d)),'FaceAlpha',1,'EdgeAlpha',0.5);
xlabel('Velocity(m/s)'); ylabel('Delay(Ts s)'); zlabel('Gain(dB)');

%% 1d cost function
Pn= findPn(Rxx_theor,M);
del_est=[140,110,30];
vel_est=[20,66,120];
[cost1d]=OneDCost(del_est,vel_est,Fjvec,r,Pn,J,c(:,1),Nsc,K);

figure;
plot(20*log10(cost1d));
xlabel('DOA(degrees)'); ylabel('Gain(dB)');
[~,DOAest]=maxk(cost1d,K);
%% gamma_kj estimation


