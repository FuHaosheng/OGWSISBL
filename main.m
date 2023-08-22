clear all;
close all;
clc;
%% Linear Arrays ULAs
M=7;   % Numbers of elements
position = [-(M-1)/2:(M-1)/2];
%%%%%%%%%%%%%%%%%% Adjustable Parameters %%%%%%%%%%%%%%%%%%%%%%
Mode = 'UncoherentSource';% 'CoherentSource' 'UncoherentSource'
SNR = 10;   resolution = 5;% grid interval
Snap = 1;% Number of snapshots
tol = 1e-4;maxiter = 500;
%% DOA Estimation
%% generate the signal
switch Mode
    case 'UncoherentSource'
        alpha1 = -35 + (10)*rand(1,1);
        alpha2 = -5 + (10)*rand(1,1); alpha3 = 20 + (10)*rand(1,1);
        TrueDOAs = round([alpha1 alpha2 alpha3]*100)/100;% alpha3
    case 'CoherentSource'
        alpha1 = -15 + (10)*rand(1,1);
        alpha2 = 10 + (10)*rand(1,1); 
        TrueDOAs = round([alpha1 alpha2]*100)/100;% alpha3
end
[Y] = signal(position, TrueDOAs, SNR, Snap, Mode);
%% %%%%%%%%%%%%%% the proposed method (2022) %%%%%%%%%%%%%%%%%%
tim0 = clock;
paras.Y = Y; paras.resolution = resolution; paras.position = position;
paras.tol = tol; paras.maxiter = maxiter;
[Pm_wsi,search_area_wsi] = OGWSISBL(paras);
tim1 = clock;
[~, Pm_wsi_id] = findpeaks((Pm_wsi),'sortstr','descend');
KP_wsi = min( length(TrueDOAs), length(Pm_wsi_id) );
Pm_wsi_ID = sort( Pm_wsi_id(1:KP_wsi), 'ascend' );
DOA_wsi_Est = search_area_wsi( Pm_wsi_ID(1:KP_wsi) );
RMSE = sqrt( ( norm( DOA_wsi_Est.' - TrueDOAs(1:KP_wsi) ).^2 ) / KP_wsi );
time = etime(tim1,tim0);
%--------------------------------------------------------------------------
%------------------------------- Figure -----------------------------------
%--------------------------------------------------------------------------
LineWidth = 1.5;
figure('NumberTitle','off','Name','DOAs Estimation')
plot(TrueDOAs,zeros(1,length(TrueDOAs)),'mo','Linewidth',LineWidth,'MarkerSize',5);
hold on; grid on; plot(search_area_wsi,10*log10(Pm_wsi/max(Pm_wsi)),'k-','Linewidth',LineWidth); hold off;
xlabel('\theta(degree)','FontName','Times','FontSize',15); ylabel('Normalized Power(dB)','FontName','Times','FontSize',15);
title( { strcat( 'SNR is', 32, num2str(SNR), 32, 'dB,', 32, ...
    'Grid interval is', 32, num2str(resolution), 32, 'deg and', 32, 'Snapshot number is', 32, num2str(Snap) ); ...
    strcat( 'RMSE is', 32, num2str(RMSE), 32, 'and', 32, 'time is', 32, num2str(time), 32, 'sec' )} );
legend('True DOAs','Proposed OGWSISBL');



