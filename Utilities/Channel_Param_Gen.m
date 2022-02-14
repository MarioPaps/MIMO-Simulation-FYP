%This function generates
%beta coefficients for all users and subcarriers
%Doppler velocities for all users and paths
%delays,DODs,DOAs  for the desired user
function[delays,beta,DODs,DOAs,vels]= Channel_Param_Gen()
    %declare constants
    K=3;
    Nsc=5;
    M=5;
    %fading coefficients
    %3 rows and 5 columns
    betas1=[0.4,0.5,0.8,0.2,0.8*exp(-1i*deg2rad(65));
            0.75, 0.32,0.65,0.88, 0.54; 
            0.67,0.75, 0.6,0.49,0.46];
    betas2= [0.23,0.48,0.76, 0.45,0.5*exp(1i*deg2rad(36));
            0.45, 0.34,0.86,0.36, 0.46;
            0.64,0.82,0.43,0.85,0.91]; 
    betas3= [0.29,0.44,0.72, 0.41,0.5*exp(1i*deg2rad(39)); 0.49, 0.33, 0.81, 0.46, 0.41;
            0.62,0.84, 0.46,0.75,0.88]; 
    betas4= [0.25,0.43,0.74, 0.49,0.51*exp(1i*deg2rad(16)); 0.44, 0.32,0.84,0.39 0.45;
            0.62,0.80, 0.49,0.81,0.72]; 
    betas5= [0.35,0.49,0.58, 0.35,0.52*exp(1i*deg2rad(46)); 0.39, 0.59,0.81,0.49, 0.31;
            0.78,0.72,0.81,0.57,0.76]; 
    
    beta=cell(1,M);
    beta{1}=betas1;
    beta{2}=betas2;
    beta{3}=betas3;
    beta{4}=betas4;
    beta{5}=betas5;

    delays= zeros(M,K);
    delays(1,:)= [140 110 30];
    temp_del=[120,45,60; 90,50,70; 30,55,80; 40,85,100];
    delays(2:end,:)=temp_del;

    DODs= zeros(M,K);
    DODs(1,:)= [35 40  60];
    DODs(2:end,:)=[87 49 67; 110 172 148; 75 83 85; 97 12 87];
   

    DOAs= zeros(M,K);
    DOAs(1,:)= [60 200 280];
    DOAs(2:end,:)=[13 76 98; 83 93 102; 115 32 83; 39 73 87];
    
    
    vels= zeros(M,K);
    vels(1,:)= [20 66 120];
    vels(2:end,:)= [31 79 94; 83 87 101; 40 31 43; 29 45 86];
end



