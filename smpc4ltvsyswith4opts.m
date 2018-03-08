% unconstrianed mpc of linear time varying system
% trajectory tracking simulation
% comparison of 5 different optimization methods
clear all; clc; close all;
pathtype = 'lorenz';   %   lorenz  duffing

switch pathtype
    case 'lorenz';
        obj1 = load('sim_lorenz_fminsearch');
        obj2 = load('sim_lorenz_fminlbfgs');
        obj3 = load('sim_lorenz_fminlbfgs_grad');
        obj4 = load('sim_lorenz_fmincon');
        obj5 = load('sim_lorenz_minimize');    
    case 'duffing';
        obj1 = load('sim_duffing_fminsearch');
        obj2 = load('sim_duffing_fminlbfgs');
        obj3 = load('sim_duffing_fminlbfgs_grad');
        obj4 = load('sim_duffing_fmincon');
        obj5 = load('sim_duffing_minimize');          
    otherwise;
end

T = obj1.T;
refpath = obj1.refpred;

linewidth = 1.2;
yyl = {'Y_1(k)', 'Y_2(k)'};
for i=1:2
    figure(i);
    plot(T, refpath(T,i), 'r:', 'linewidth', linewidth); hold on;
    plot(obj1.gppred(:,i), 'b-', 'linewidth', linewidth); hold on;
    plot(obj2.gppred(:,i), 'k-', 'linewidth', linewidth); hold on;
    plot(obj3.gppred(:,i), 'b.-', 'linewidth', linewidth); hold on;
    plot(obj4.gppred(:,i), 'k.-', 'linewidth', linewidth); hold on;
    plot(obj5.gppred(:,i), 'k--', 'linewidth', linewidth);
    ylabel(yyl{i});
    xlim([0, numel(T)]);
    xlabel('sampling time k');
    legend('Reference', 'fminsearch', 'fminlbfgs',...
        'fminlbfgs-grad', 'fmincon', 'minimize');
end


figure(3);
plot(obj1.fval, 'k-', 'linewidth', linewidth); hold on;
plot(obj2.fval, 'b-', 'linewidth', linewidth); hold on;
plot(obj3.fval, 'r-', 'linewidth', linewidth); hold on;
plot(obj4.fval, 'k+-', 'linewidth', linewidth); hold on;
plot(obj5.fval, 'r.-', 'linewidth', linewidth);
xlim([0, numel(T)]);
xlabel('sampling time k');
ylabel('Cost Values');
legend('fminsearch', 'fminlbfgs',...
    'fminlbfgs-grad', 'fmincon', 'minimize');

figure(4);
plot(obj1.fval, 'k-', 'linewidth', linewidth); hold on;
plot(obj4.fval, 'b-', 'linewidth', linewidth);
xlim([0, numel(T)]);
xlabel('sampling time k');
ylabel('Cost Values');
legend('fminsearch', 'minimize');

figure(5);
plot(obj1.fval-obj4.fval, 'b-', 'linewidth', linewidth);
xlim([0, numel(T)]);
xlabel('sampling time k');
ylabel('Cost Values');
legend('fminsearch-minimize');


