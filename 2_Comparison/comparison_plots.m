%code to plot the mat files obtaining from running the comparison code
%receiver comparison plot
load("doppstar_result.mat");
load("rake_resultfinal.mat");
load("spacesub_result2.mat");

NFR= 10.^((0:5:60)./10); %NFR levels
figure;
plot(10*log10(NFR),mean(doppstar_result),'-o','DisplayName','Doppler STAR Subspace Rx','LineWidth',2);
hold on;
plot(10*log10(NFR),mean(spacesub_result),'-o','Color',[107 76 154]./255,'DisplayName','Massive MIMO Subspace Rx','LineWidth',2);
plot(10*log10(NFR),(rake_result(1,:)),'-o','Color','r','DisplayName','Doppler STAR RAKE Rx','LineWidth',2);
xlabel('NFR (dB)'); ylabel('SNIRout_{1} (dB)');
title('SNIRout Comparison');grid on; grid minor;
yticks(-50:10:50);
ax = gca; 
ax.FontSize = 11; 
legend('show');
%% subcarriers plot
load("SNIR_sub.mat");
figure;
plot((1:width(SNIR_sub)),mean(10*log10(real(SNIR_sub))),'-o','LineWidth',1.5);
grid on; grid minor;
xlabel('Number of Subcarriers'); ylabel('SNIRout_{1} (dB)'); 
title('Impact of {\it N_{sc}} Variation on SNIRout_{j} ');
ax = gca; 
ax.FontSize = 11; 