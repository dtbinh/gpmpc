clear;clc;close all;
set(0, 'DefaultAxesFontName', 'Arial');  
set(0, 'DefaultTextFontName', 'Arial');  


tempds1_1 = load('lorenz_H10_GPMPC1_60.mat');
tempds1_2 = load('lorenz_H10_GPMPC1_80.mat');
tempds1_3 = load('lorenz_H10_GPMPC1_100.mat');

tempds2_1 = load('lorenz_H10_GPMPC2_60.mat');
tempds2_2 = load('lorenz_H10_GPMPC2_80.mat');
tempds2_3 = load('lorenz_H10_GPMPC2_100.mat');

linewith = 1;
for i=1:2
    figure(i);
    plot(tempds1_1.refpred(:, 2*i-1), 'k:', 'linewidth', 1.5); hold on;
    plot(tempds1_1.gppred(:, 2*i-1), 'r--', 'linewidth', linewith); hold on;
    plot(tempds1_2.gppred(:, 2*i-1), '--', 'linewidth', linewith,...
        'color', [0 100 0]/255); hold on;
    plot(tempds1_3.gppred(:, 2*i-1), 'k-', 'linewidth', linewith); hold on;
    plot(tempds2_1.gppred(:, 2*i-1), 'r:', 'linewidth', linewith); hold on;
    plot(tempds2_2.gppred(:, 2*i-1), ':', 'linewidth', linewith,...
        'color', [0 100 0]/255); hold on;
    plot(tempds2_3.gppred(:, 2*i-1), 'b-', 'linewidth', linewith); hold on;
    xlim([0 180]);
    xlabel('sample time k');
    ylabel(strcat('Y_',int2str(i),'(k)'));
    legend('Reference', 'GPMPC1-60%', 'GPMPC1-80%','GPMPC1-90%',...
        'GPMPC2-60%','GPMPC2-80%','GPMPC2-90%');
end


