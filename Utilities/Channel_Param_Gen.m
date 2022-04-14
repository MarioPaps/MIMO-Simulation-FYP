%This function generates
%beta coefficients for all users and subcarriers
%Doppler velocities for all users and paths
%delays,DODs,DOAs  for the desired user
function[delays,beta,DODs,DOAs,vels]= Channel_Param_Gen(set_same_vel,set_same_DOA)
    %declare constants
    K=3;
    M=5;
    %fading coefficients
    %3 rows and 5 columns
    %j=1, j=2... of user 1
%     betas1=[0.75*exp(1i*deg2rad(78)) 0.66 0.5; 0.65*exp(1i*deg2rad(32)) 0.48 0.61; 0.62*exp(1i*deg2rad(64)) 0.37 0.41; 0.54*exp(1i*deg2rad(-40)) 0.38 0.39; 0.27*exp(1i*deg2rad(-75)) 0.44 0.28];
%     betas1=betas1';
    %betas for regular plots
    betas1=[0.82*exp(1i*deg2rad(78)) 0.75 0.67; 0.9*exp(1i*deg2rad(32)) 0.88 0.61; 0.82*exp(1i*deg2rad(64)) 0.67 0.71; 0.84*exp(1i*deg2rad(-40)) 0.88 0.74; 0.87*exp(1i*deg2rad(-75)) 0.84 0.78];
    betas1=betas1.';
    betas2=[0.1 0.11 0.12; 0.23 0.01 0.06; 0.08 0.10 0.15; 0.12 0.23 0.26; 0.27 0.23 0.28];
    betas2=betas2.';
    betas3=[0.03 0.08 0.04; 0.12 0.14 0.04; 0.17 0.24 0.29; 0.18 0.22 0.27; 0.21 0.24 0.29];
    betas3=betas3.';
    betas4=[0.02 0.21 0.07; 0.14 0.18 0.21; 0.25 0.22 0.19; 0.20 0.24 0.31; 0.32 0.35 0.39];
    betas4=betas4.';
    betas5=[0.25 0.18 0.14; 0.19 0.24 0.03; 0.05 0.03 0.10; 0.12 0.14 0.11; 0.17 0.12 0.15];
    betas5=betas5.';

    %betas that worked 27 Feb
%      betas1=[0.82*exp(1i*deg2rad(78)) 0.75 0.67; 0.9*exp(1i*deg2rad(32)) 0.88 0.61; 0.82*exp(1i*deg2rad(64)) 0.67 0.71; 0.84*exp(1i*deg2rad(-40)) 0.88 0.74; 0.87*exp(1i*deg2rad(-75)) 0.84 0.78];
%     betas1=betas1';

%     %betas for Monte Carlo simulation
%     betas1=[0.43*exp(1i*deg2rad(78)) 0.15 0.37; 0.42*exp(1i*deg2rad(32)) 0.49 0.5; 0.42*exp(1i*deg2rad(64)) 0.37 0.48; 0.48*exp(1i*deg2rad(-40)) 0.29 0.36; 0.36*exp(1i*deg2rad(-75)) 0.35 0.26];
%     betas1=betas1.';
%     betas2=[0.1 0.11 0.12; 0.23 0.01 0.06; 0.08 0.10 0.15; 0.12 0.23 0.26; 0.27 0.23 0.28];
%     betas2=betas2.';
%     betas3=[0.03 0.08 0.04; 0.12 0.14 0.04; 0.17 0.24 0.29; 0.18 0.22 0.27; 0.21 0.24 0.29];
%     betas3=betas3.';
%     betas4=[0.02 0.21 0.07; 0.14 0.18 0.21; 0.25 0.22 0.19; 0.20 0.24 0.31; 0.32 0.35 0.39];
%     betas4=betas4.';
%     betas5=[0.25 0.18 0.14; 0.19 0.24 0.03; 0.05 0.03 0.10; 0.12 0.14 0.11; 0.17 0.12 0.15];
%     betas5=betas5.';





%      betas1=[0.82,0.76,0.8,0.2,0.8*exp(-1i*deg2rad(65));
%             0.75, 0.32,0.65,0.88, 0.54; 
%             0.67,0.75, 0.6,0.49,0.46];
%     betas2= [0.23,0.48,0.76, 0.45,0.5*exp(1i*deg2rad(36));
%             0.45, 0.34,0.86,0.36, 0.46;
%             0.64,0.82,0.43,0.85,0.91]; 
%     betas3= [0.29,0.44,0.72, 0.41,0.5*exp(1i*deg2rad(39)); 0.49, 0.33, 0.81, 0.46, 0.41;
%             0.62,0.84, 0.46,0.75,0.88]; 
%     betas4= [0.25,0.43,0.74, 0.49,0.51*exp(1i*deg2rad(16)); 0.44, 0.32,0.84,0.39 0.45;
%             0.62,0.80, 0.49,0.81,0.72]; 
%     betas5= [0.35,0.49,0.58, 0.35,0.52*exp(1i*deg2rad(46)); 0.39, 0.59,0.81,0.49, 0.31;
%             0.78,0.72,0.81,0.57,0.76]; 
%     betas1=[0.74,0.95,0.8,0.87,0.8*exp(-1i*deg2rad(65));
%             0.75, 0.92,0.65,0.88, 0.59; 
%             0.67,0.75, 0.69,0.92,0.63];
%     betas2= [0.23,0.01,0.16, 0.04,0.05*exp(1i*deg2rad(36));
%             0.08, 0.14,0.06,0.07, 0.04;
%             0.04,0.02,0.03,0.05,0.01]; 
%     betas3= [0.29,0.02,0.01, 0.14,0.07*exp(1i*deg2rad(39)); 0.02, 0.03, 0.13, 0.16, 0.05;
%             0.062,0.041, 0.016,0.015,0.08]; 
%     betas4= [0.25,0.43,0.02, 0.09,0.15*exp(1i*deg2rad(16)); 0.14, 0.22,0.09,0.03 0.01;
%             0.03,0.10, 0.19,0.11,0.12]; 
%     betas5= [0.001,0.049,0.008, 0.05,0.25*exp(1i*deg2rad(46)); 0.01, 0.02,0.03,0.05, 0.12;
%             0.18,0.12,0.03,0.06,0.11]; 
    
%     beta=cell(1,M);
%     beta{1}=betas1;
%     beta{2}=betas2;
%     beta{3}=betas3;
%     beta{4}=betas4;
%     beta{5}=betas5;

    beta=[betas1, betas2,betas3,betas4,betas5];

    delays= zeros(M,K);
    delays(1,:)= [140 110 30];
    temp_del=[120,45,60; 90,50,70; 30,55,80; 40,85,100];
    delays(2:end,:)=temp_del;

    DODs= zeros(M,K);
    DODs(1,:)= [35 40  60];
    DODs(2:end,:)=[87 18 137; 110 172 148; 75 83 85; 97 12 81];
   

    DOAs= zeros(M,K);
    DOAs(1,:)= [60 200 280];
    DOAs(2:end,:)=[13 76 98; 83 93 102; 115 32 83; 39 73 87];
    
    
    vels= zeros(M,K);
    vels(1,:)= [20 66 120];
    vels(2:end,:)= [31 79 94; 83 87 101; 40 31 43; 29 45 86];
    
    if(set_same_vel==1)
       vels(1,:)=20;
    end
    if(set_same_DOA==1)
        DOAs(1,2)=DOAs(1,3); %make 2 paths have DOA 280
        DOAs(2,2)= DOAs(1,3); %user 2 on path 2 also has DOA 280
    end
end



