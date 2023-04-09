%channel parameters function
%set_same_vel,set_same_DOA,Nsc_comp,Nsc
function [delays,beta,DODs,DOAs,vels] = ChannelParam(set_same_vel,set_same_DOA,Nsc_comp,Nsc)
    M=5;
    K=3;
    %fading coefficients- max size - K*(M*Nsc) ->max values:K=3, M=5, Nsc=10
    %user 1 betas
  %  betas1=0.6*ones(K,M*Nsc);
   betas1=[0.52*exp(-1i*deg2rad(78)) 0.75 0.67 0.49 0.58 0.78 0.74 0.44 0.5 0.63;
         0.6*exp(-1i*deg2rad(32)) 0.88 0.61 0.42 0.39 0.43 0.53 0.61 0.65 0.41;
         0.42*exp(-1i*deg2rad(64)) 0.67 0.71 0.49*exp(-1i*deg2rad(-40)) 0.58 0.74 0.55 0.73 0.65 0.69];
    %user 2 betas
    betas2=[0.08 0.01 0.05 0.03 0.01 0.06 0.08 0.10 0.15 0.03 ;
        0.26 0.29 0.15*exp(-1i*deg2rad(12)) 0.12 0.14 0.09 0.26 0.24 0.28 0.15;
        0.16 0.19 0.08 0.16 0.19 0.23 0.2 0.18 0.13 0.19];
    %user 3 betas
    betas3=[0.15 0.13 0.19 0.28 0.03 0.06 0.03 0.12 0.13 0.24 ;
    0.21 0.10 0.19 0.3 0.29 0.13 0.19 0.12 0.24 0.16;
    0.16 0.1 0.14 0.19 0.18 0.24 0.18 0.28 0.26 0.25];
    %user 4 betas
    betas4=[0.25 0.12 0.11 0.29 0.12 0.03 0.09 0.15 0.23 0.19 ;
    0.15 0.19 0.21 0.3 0.31 0.18 0.12 0.14 0.17 0.23;
    0.11 0.18 0.15 0.12 0.2 0.14 0.21 0.1 0.18 0.25];
    %user 5 betas
    betas5=[0.06 0.22 0.15 0.14 0.32 0.13 0.19 0.15 0.18 0.29 ;
    0.21 0.22 0.21 0.3 0.13 0.22 0.25 0.26 0.27 0.12;
    0.18 0.21 0.24 0.22 0.21 0.18 0.15 0.09 0.22 0.1];
    %assign
    beta=[betas1,betas2,betas3,betas4,betas5];
    if(Nsc<10)
        betas1= betas1(:,1:Nsc);
        betas2= betas2(:,1:Nsc);
        betas3= betas3(:,1:Nsc);
        betas4= betas4(:,1:Nsc);
        betas5= betas5(:,1:Nsc);
        beta=[betas1,betas2,betas3,betas4,betas5];

       
    end
    
    %delays
    delays= zeros(M,K);
    delays(1,:)= [140 110 30];
    delays(2:end,:)=[20,45,60; 90,50,70; 30,55,80; 40,85,100];
    if(Nsc<5) || (Nsc_comp==1)
        delays=floor(delays/10);
    end

    %Doppler velocities
    vels= zeros(M,K);
    vels(1,:)= [20 66 120];
    vels(2:end,:)= [31 79 94; 83 87 101; 40 31 43; 29 45 86];
    if(set_same_vel==1)
       vels(1,:)=20;
    end
    
    %DODs
    DODs= zeros(M,K);
    DODs(1,:)= [35 40  28];
    DODs(2:end,:)=[87 18 137; 110 172 148; 75 83 85; 97 12 81];

    %DOAs
    DOAs= zeros(M,K);
    DOAs(1,:)= [60 200 280];
    DOAs(2:end,:)=[13 76 98; 83 93 102; 115 32 83; 39 73 87];
    if(set_same_DOA==1)
        DOAs(1,2)=DOAs(1,3); %make 2 paths of user 1 have DOA 280
        DOAs(2,2)= DOAs(1,3); %user 2 on path 2 also has DOA 280
    end

end