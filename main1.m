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
paras.tol = tol; paras.maxiter = maxiter;  paras.etc = M; paras.r_step = 0.001;
[Pm_wsi,search_area_wsi] = OGWSISBL(paras);
tim1 = clock;
[~, Pm_wsi_id] = findpeaks((Pm_wsi),'sortstr','descend');
KP_wsi = min( length(TrueDOAs), length(Pm_wsi_id) );
Pm_wsi_ID = sort( Pm_wsi_id(1:KP_wsi), 'ascend' );
DOA_wsi_Est = search_area_wsi( Pm_wsi_ID(1:KP_wsi) );
RMSE = sqrt( ( norm( DOA_wsi_Est.' - TrueDOAs(1:KP_wsi) ).^2 ) / KP_wsi );
time = etime(tim1,tim0);
%% OGSBI (2013)
tim2 = clock;
[Pm_ogsbi,search_area_ogsbi]=OGSBI(paras);
tim3 = clock;
[~, Pm_ogsbi_id] = findpeaks((Pm_ogsbi),'sortstr','descend');
KP_ogsbi = min(length(TrueDOAs),length(Pm_ogsbi_id));
Pm_ogsbi_ID = Pm_ogsbi_id(1:KP_ogsbi);
Pm_ogsbi_ID = sort(Pm_ogsbi_ID,'ascend');
DOA_ogsbi_Est = search_area_ogsbi(Pm_ogsbi_ID(1:KP_ogsbi));
RMSE_OGSBI = sqrt( (norm(DOA_ogsbi_Est - TrueDOAs(1:KP_ogsbi)).^2)/KP_ogsbi );
time_ogsbi = etime(tim3,tim2);
%% PSBL (2016)
tim4 = clock;
[Pm_psbl,search_area_psbl]=PSBL(paras);
tim5 = clock;
[~, Pm_psbl_id] = findpeaks((Pm_psbl),'sortstr','descend');
KP_psbl = min(length(TrueDOAs),length(Pm_psbl_id));
Pm_psbl_ID = Pm_psbl_id(1:KP_psbl);
Pm_psbl_ID = sort(Pm_psbl_ID,'ascend');
DOA_psbl_Est = search_area_psbl(Pm_psbl_ID(1:KP_psbl));
RMSE_PSBL = sqrt( (norm(DOA_psbl_Est - TrueDOAs(1:KP_psbl)).^2)/KP_psbl );
time_psbl = etime(tim5,tim4);
%% RSBL (2017)
tim6 = clock;
[Pm_rsbl,search_area_rsbl] = RSBL(paras);
tim7 = clock;
[~, Pm_rsbl_id] = findpeaks((Pm_rsbl),'sortstr','descend');
KP_rsbl = min(length(TrueDOAs),length(Pm_rsbl_id));
Pm_rsbl_ID = Pm_rsbl_id(1:KP_rsbl);
Pm_rsbl_ID = sort(Pm_rsbl_ID,'ascend');
DOA_rsbl_Est = search_area_rsbl(Pm_rsbl_ID(1:KP_rsbl));
RMSE_RSBL = sqrt( (norm(DOA_rsbl_Est - TrueDOAs(1:KP_rsbl)).^2)/KP_rsbl );
time_rsbl = etime(tim7,tim6);
%% Real-VBI (2021)
tim8 = clock;
[Pm_rvbi,search_area_rvbi]=real_VBI(paras); 
tim9 = clock;
[~, Pm_rvbi_id] = findpeaks((Pm_rvbi),'sortstr','descend');
KP_rvbi = min(length(TrueDOAs),length(Pm_rvbi_id));
Pm_rvbi_ID = Pm_rvbi_id(1:KP_rvbi);
Pm_rvbi_ID = sort(Pm_rvbi_ID,'ascend');
DOA_rvbi_Est = search_area_rvbi(Pm_rvbi_ID(1:KP_rvbi));
RMSE_RVBI = sqrt( (norm(DOA_rvbi_Est - TrueDOAs(1:KP_rvbi)).^2)/KP_rvbi );
time_rvbi = etime(tim9,tim8);
%% CGDP (2021)
tim10 = clock;
[Pm_cgdp1,search_area_cgdp1]=CGDP_SBL1(paras);
tim11 = clock;
[~, Pm_cgdp_id1] = findpeaks((Pm_cgdp1),'sortstr','descend');
KP_cgdp1 = min(length(TrueDOAs),length(Pm_cgdp_id1));
Pm_cgdp_ID1 = Pm_cgdp_id1(1:KP_cgdp1);
Pm_cgdp_ID1 = sort(Pm_cgdp_ID1,'ascend');
DOA_cgdp_Est1 = search_area_cgdp1(Pm_cgdp_ID1(1:KP_cgdp1));
RMSE_CGDP1 = sqrt( (norm(DOA_cgdp_Est1 - TrueDOAs(1:KP_cgdp1)).^2)/KP_cgdp1 );
time_cgdp1 = etime(tim11,tim10);
%%
%--------------------------------------------------------------------------
%------------------------------- Figure -----------------------------------
%--------------------------------------------------------------------------
LineWidth = 1.5;
figure('NumberTitle','off','Name','DOAs Estimation')
subplot(3,2,1);plot(TrueDOAs,zeros(1,length(TrueDOAs)),'mo','Linewidth',LineWidth,'MarkerSize',5);
hold on; grid on; plot(search_area_wsi,10*log10(Pm_wsi/max(Pm_wsi)),'k-','Linewidth',LineWidth); hold off;
xlabel('\theta(degree)','FontName','Times','FontSize',15); ylabel('Normalized Power(dB)','FontName','Times','FontSize',15);
title( { strcat( 'SNR is', 32, num2str(SNR), 32, 'dB,', 32, ...
    'Grid interval is', 32, num2str(resolution), 32, 'deg and', 32, 'Snapshot number is', 32, num2str(Snap) ); ...
    strcat( 'RMSE is', 32, num2str(RMSE), 32, 'and', 32, 'time is', 32, num2str(time), 32, 'sec' )} );
legend('True DOAs','Proposed OGWSISBL');

subplot(3,2,2);plot(TrueDOAs,zeros(1,length(TrueDOAs)),'mo','Linewidth',LineWidth,'MarkerSize',5);
hold on; grid on; plot(search_area_ogsbi,10*log10(Pm_ogsbi/max(Pm_ogsbi)),'r-','Linewidth',LineWidth); hold off;
xlabel('\theta(degree)','FontName','Times','FontSize',15); ylabel('Normalized Power(dB)','FontName','Times','FontSize',15);
title( { strcat( 'SNR is', 32, num2str(SNR), 32, 'dB,', 32, ...
    'Grid interval is', 32, num2str(resolution), 32, 'deg and', 32, 'Snapshot number is', 32, num2str(Snap) ); ...
    strcat( 'RMSE is', 32, num2str(RMSE_OGSBI), 32, 'and', 32, 'time is', 32, num2str(time_ogsbi), 32, 'sec' )} );
legend('True DOAs','OGSBI');

subplot(3,2,3);plot(TrueDOAs,zeros(1,length(TrueDOAs)),'mo','Linewidth',LineWidth,'MarkerSize',5);
hold on; grid on; plot(search_area_psbl,10*log10(Pm_psbl/max(Pm_psbl)),'k-','Linewidth',LineWidth); hold off;
xlabel('\theta(degree)','FontName','Times','FontSize',15); ylabel('Normalized Power(dB)','FontName','Times','FontSize',15);
title( { strcat( 'SNR is', 32, num2str(SNR), 32, 'dB,', 32, ...
    'Grid interval is', 32, num2str(resolution), 32, 'deg and', 32, 'Snapshot number is', 32, num2str(Snap) ); ...
    strcat( 'RMSE is', 32, num2str(RMSE_PSBL), 32, 'and', 32, 'time is', 32, num2str(time_psbl), 32, 'sec' )} );
legend('True DOAs','PSBL');

subplot(3,2,4);plot(TrueDOAs,zeros(1,length(TrueDOAs)),'mo','Linewidth',LineWidth,'MarkerSize',5);
hold on; grid on; plot(search_area_rsbl,10*log10(Pm_rsbl/max(Pm_rsbl)),'k-','Linewidth',LineWidth); hold off;
xlabel('\theta(degree)','FontName','Times','FontSize',15); ylabel('Normalized Power(dB)','FontName','Times','FontSize',15);
title( { strcat( 'SNR is', 32, num2str(SNR), 32, 'dB,', 32, ...
    'Grid interval is', 32, num2str(resolution), 32, 'deg and', 32, 'Snapshot number is', 32, num2str(Snap) ); ...
    strcat( 'RMSE is', 32, num2str(RMSE_RSBL), 32, 'and', 32, 'time is', 32, num2str(time_rsbl), 32, 'sec' )} );
legend('True DOAs','RSBL');

subplot(3,2,5);plot(TrueDOAs,zeros(1,length(TrueDOAs)),'mo','Linewidth',LineWidth,'MarkerSize',5);
hold on; grid on; plot(search_area_rvbi,10*log10(Pm_rvbi/max(Pm_rvbi)),'k-','Linewidth',LineWidth); hold off;
xlabel('\theta(degree)','FontName','Times','FontSize',15); ylabel('Normalized Power(dB)','FontName','Times','FontSize',15);
title( { strcat( 'SNR is', 32, num2str(SNR), 32, 'dB,', 32, ...
    'Grid interval is', 32, num2str(resolution), 32, 'deg and', 32, 'Snapshot number is', 32, num2str(Snap) ); ...
    strcat( 'RMSE is', 32, num2str(RMSE_RVBI), 32, 'and', 32, 'time is', 32, num2str(time_rvbi), 32, 'sec' )} );
legend('True DOAs','Real Value-BI');

subplot(3,2,6);plot(TrueDOAs,zeros(1,length(TrueDOAs)),'mo','Linewidth',LineWidth,'MarkerSize',5);
hold on; grid on; plot(search_area_cgdp1,10*log10(Pm_cgdp1/max(Pm_cgdp1)),'k-','Linewidth',LineWidth); hold off;
xlabel('\theta(degree)','FontName','Times','FontSize',15); ylabel('Normalized Power(dB)','FontName','Times','FontSize',15);
title( { strcat( 'SNR is', 32, num2str(SNR), 32, 'dB,', 32, ...
    'Grid interval is', 32, num2str(resolution), 32, 'deg and', 32, 'Snapshot number is', 32, num2str(Snap) ); ...
    strcat( 'RMSE is', 32, num2str(RMSE_CGDP1), 32, 'and', 32, 'time is', 32, num2str(time_cgdp1), 32, 'sec' )} );
legend('True DOAs','CGDP-SBL1');

