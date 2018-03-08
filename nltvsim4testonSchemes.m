clear;clc;close all;
set(0, 'DefaultAxesFontName', 'Arial');  
set(0, 'DefaultTextFontName', 'Arial');  

tem = 'lorenz';     %  lorenz    step

switch tem
    case 'lorenz';
        tempds1 = load('lorenz_H1_fmincon.mat');
        tempds2 = load('lorenz_H1_gpmpc1.mat');
        tempds3 = load('lorenz_H1_gpmpc2.mat');
    case 'step';
        tempds1 = load('step_H1_fmincon.mat');
        tempds2 = load('step_H1_gpmpc1.mat');
        tempds3 = load('step_H1_gpmpc2.mat');
end

linewith = 1;
for i=1:2
    figure(i);
    plot(tempds1.refpred(:, 2*i-1), 'k:', 'linewidth', 1.5); hold on;
    plot(tempds1.gppred(:, 2*i-1), 'r--', 'linewidth', linewith); hold on;
    plot(tempds2.gppred(:, 2*i-1), 'k-', 'linewidth', linewith); hold on;
    plot(tempds3.gppred(:, 2*i-1), 'b-', 'linewidth', linewith); hold on;
    xlim([0 189]);
    xlabel('sample time k');
    ylabel(strcat('Y_',int2str(i),'(k)'));
    legend('Reference', 'Nonlinear GPMPC', 'GPMPC1','GPMPC2');
end