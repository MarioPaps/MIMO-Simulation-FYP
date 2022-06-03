%code to plot the mat files obtaining from running the comparison code
%receiver comparison plot
NFR= 10.^((0:5:60)./10); %NFR levels
figure;
plot(10*log10(NFR),mean(doppstar_result),'-o','DisplayName','Doppler STAR Subspace Rx');
hold on;
plot(10*log10(NFR),(rake_result(1,:)),'-o','DisplayName','Doppler STAR RAKE Rx');
plot(10*log10(NFR),mean(spacesub_result),'-o','Color','g','DisplayName','maMIMO Subspace Rx');
xlabel('NFR(dB)'); ylabel('SNIRout_{1} (dB)');
title('SNIRout Comparison');grid on; grid minor;
yticks(-50:10:50);
ax = gca; 
ax.FontSize = 11; 
legend('show');
%% subcarriers plot
figure;
plot((1:width(SNIR_sub)),mean(10*log10(real(SNIR_sub))),'-o');
grid on; grid minor;
xlabel('Number of Subcarriers'); ylabel('SNIRout_{1} (dB)'); 
title('Impact of N_{sc} Variation on SNIRout_{j}');
ax = gca; 
ax.FontSize = 11; 