% SNR=[0.5,3.1623,5,40,50,250,500,25000,5000];
% L=200;
% xaxis= L*SNR;
% load("rmse.mat");
% mean_rmse=mean(RMSEradvel(1:82,1:9));
% temp=mean_rmse(9);
% mean_rmse(9)= mean_rmse(8);
% mean_rmse(8)= temp;
% % mean_rmse(8)= mean_rmse(8)/10;
% figure;
% loglog(xaxis(1:9),mean_rmse,'-s','DisplayName','Radial Vel RMSE');
% % line(xaxis,RMSEradvel);
% xlabel('SNRxL'); ylabel('Estimation RMSE');
% grid on;
% hold off;
% ylim([10e-5 10e1]); legend('show');
addpath('.\1_MC Results\')

array1=load("rmse.mat");
array2=load("newest.mat");
array3=load("greatestall.mat");

array1=cell2mat(struct2cell(array1));
array2=cell2mat(struct2cell(array2));
array3=cell2mat(struct2cell(array3));

RMSEradvel= zeros(100,9);
RMSEradvel(1:82,1:7)= array1(1:82,1:7);
RMSEradvel(:,8:9)= array2(:,8:9);
RMSEradvel(83:end,:)=array3(83:end,:);

SNR=[0.5,3.1623,5,40,50,250,500,2500,5000];
L=200;
xaxis= L*SNR;
mean_rmse=mean(RMSEradvel(1:100,1:9));

figure;
loglog(xaxis(1:9),mean_rmse,'-o','DisplayName','Radial Velocity (m/s) RMSE');
% line(xaxis,RMSEradvel);
xlabel('SNRxL'); ylabel('Estimation RMSE');
grid on;
hold on;
loglog(xaxis(1:9),ress,'-s','DisplayName','DOA RMSE');
set(gca,'YLim',[10e-7 10e1],'YTick',10.^(-7:1))
% ylim([10e-6 10e1]);
legend('show');
%% velocity data
velest1=load("vel_est.mat");
velest89=load("allvels.mat");
velest1last= load("vels4ever.mat");

velest1=cell2mat(struct2cell(velest1));
velest89=cell2mat(struct2cell(velest89));
velest1last=cell2mat(struct2cell(velest1last));


vel_est=zeros(9,300);
vel_est(1:7,1:246)= velest1(1:7,1:246);
vel_est(1:7,247:300)= velest1last(1:7,:);

vel_est(8:9,:)= velest89(8:9,:); %complete summary of estimated velocities



