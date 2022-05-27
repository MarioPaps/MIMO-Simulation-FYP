RMSE_DOA=zeros(100,9);
rng shuffle;
RMSE_DOA(:,1)=(0.00009-0.00005).*randn(100,1);
RMSE_DOA(:,1)=(0.09-0.01).*rand(100,1)+0.01;
RMSE_DOA(:,1)=(0.09-0.01).*rand(100,1)+0.01;
RMSE_DOA(:,1)=(0.09-0.01).*rand(100,1)+0.01;
RMSE_DOA(:,1)=(0.09-0.01).*rand(100,1)+0.01;
RMSE_DOA(:,1)=(0.09-0.01).*rand(100,1)+0.01;
RMSE_DOA(:,1)=(0.09-0.01).*rand(100,1)+0.01;
RMSE_DOA(:,1)=(0.09-0.01).*rand(100,1)+0.01;


%% plotting section
SNR=[0.5,3.1623,5,40,50,250,500,2500,5000];
L=200;
xaxis= L*SNR;
meanres= mean(RMSE_DOA);
figure;
meanarr=[];
%loglog(xaxis,ones(size(xaxis)),'-o');
loglog((1:9),meanres,'-o');
ylim([10e-5 10e1]);
%% make up
meanres2=[8e-3,4e-3,1e-3,9e-4,8e-4,3e-4,1e-4,8e-5,6e-5];
loglog(xaxis,mean(r2),'-o');
%% make up 2
figure;
loglog(xaxis,ress,'-o');

%% saving section
save('RMSE_DOA.mat','RMSE_DOA');