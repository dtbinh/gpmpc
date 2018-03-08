function nltvsims4Journalpaper2
clear;clc;close all;
set(0, 'DefaultAxesFontName', 'Arial');  
set(0, 'DefaultTextFontName', 'Arial');  

pathtype = 'step';     % duffing  lorenz  step
T=1:189;
switch pathtype
    case 'duffing';
        temp = load('nmpc_nltv_duffing.mat'); obj1 = temp.obj;
        temp = load('gpmpc_duffing_sqp.mat'); obj2 = wrapObj(temp);
        temp = load('gpmpc_duffing_quadprog.mat'); obj3 = wrapObj(temp);
        titletex1 = '"Duffing" Trajectory Tracking';
        titletex2 = 'Control inputs in "Duffing" Trajectory Tracking';
    case 'lorenz';
        temp = load('nmpc_nltv_lorenz.mat'); obj1 = temp.obj;
        temp = load('gpmpc_lorenz_sqp.mat'); obj2 = wrapObj(temp);
        temp = load('gpmpc_lorenz_quadprog.mat'); obj3= wrapObj(temp);
        titletex1 = '"Lorenz" Trajectory Tracking';
        titletex2 = 'Control inputs in "Lorenz" Trajectory Tracking';
    case 'step';
        ds=load('ds_nmpc_nltv_step.txt');
        obj1.path =  TracjectorGenerator(1:200,4);
        obj1.x = ds(:,1:4);
        obj1.u = ds(:,5:6);
        temp = load('gpmpc_step_sqp.mat'); obj2 = wrapObj(temp);
        temp = load('gpmpc_step_quadprog.mat'); obj3 = wrapObj(temp);
        titletex1 = '"Step" Trajectory Tracking';
        titletex2 = 'Control inputs in "Step" Trajectory Tracking';
end

linewith = 1;
xlb = {'Y_1(k)', 'Y_2(k)'};
figure(1);
for i=1:2
    subplot(2,1,i);
    plot(obj1.path(:,i),'r:'); hold on;
    plot(obj1.x(:,2*(i-1)+1), '-.', 'linewidth', linewith, 'color',[96 96 96]/255);  hold on;
    plot(obj2.gppred(:,2*(i-1)+1), 'k-', 'linewidth', linewith); hold on;
    plot(obj3.gppred(:,2*(i-1)+1), 'b--', 'linewidth', linewith); hold on;
%     plot(sqpobj3.gppred(:,2*(i-1)+1), '-', 'linewidth', linewith, 'Color', [0 0 0]/255); 
    xlim([0 max(T)]);
    ylabel(xlb{i});
    if i==1
        title(titletex1);
    end
end
xlabel('sample time k');
legend('Reference','NMPC', 'GPMPC1', 'GPMPC2');

figure(2);
xlb={'u_1(k)', 'u_2(k)'};
for i=1:2
    subplot(2,1,i);
    plot(obj1.u(:,i), 'r:', 'linewidth', linewith);  hold on;
    plot(obj2.u(i,:), 'k-', 'linewidth', linewith); hold on;
    plot(obj3.u(i,:), 'b--'); 
    xlim([0 max(T)]);
    ylabel(xlb{i});
    if i==1
        title(titletex2);
    end
end
xlabel('sample time k');
legend('NMPC', 'GPMPC1', 'GPMPC2');

iae1 = abs(obj2.gppred(:,[1 3])-obj1.path(1:189,:));
iae2 = abs(obj3.gppred(:,[1 3])-obj1.path(1:189,:));
xlb = {'Y_1(k)', 'Y_2(k)'};
figure(3);
subplot(2,1,1);
plot(cumsum(iae1(:,1)), 'k--', 'linewidth', linewith); hold on;
plot(cumsum(iae2(:,1)), 'b-'); 
xlim([0 max(T)]);
ylabel(xlb{1});
title({'Integral Absolute Error (IAE)', titletex1});
subplot(2,1,2);
plot(cumsum(iae1(:,2)), 'k--', 'linewidth', linewith); hold on;
plot(cumsum(iae2(:,2)), 'b-'); 
xlim([0 max(T)]);
ylabel(xlb{2});
xlabel('sample time k');
legend('GPMPC1', 'GPMPC2');

figure(4);
subplot(2,1,1);
plot(cumsum(iae1(:,1))-cumsum(iae2(:,1)), 'k-', 'linewidth', linewith);
xlim([0 max(T)]);
ylabel(xlb{1});
title({'IAE Differences of GPMPC1 and GPMPC2',titletex1});
subplot(2,1,2);
plot(cumsum(iae1(:,2))-cumsum(iae2(:,2)), 'k-', 'linewidth', linewith);
xlim([0 max(T)]);
ylabel(xlb{1});
xlabel('sample time k');
legend('GPMPC1-GPMPC2');

% figure(4);
% subplot(2,1,1);
% plot(obj2.fval, 'k-'); hold on;
% plot(obj3.fval, 'b-.');
% subplot(2,1,2);
% plot(obj3.fval-obj2.fval, 'r--');
% myzoom;
% figure(2);
% for i=1:2
%     subplot(2,1,i);
%     plot(obj1.path(:,i),':', 'linewidth', linewith,'Color', [255 0 0]/255); hold on;
% %     plot(obj1.x(:,2*(i-1)+1), 'b-.', 'linewidth', linewith);  hold on;
%     plot(lgpobj1.gppred(:,2*(i-1)+1), 'b-', 'linewidth', linewith); hold on;
%     plot(lgpobj2.gppred(:,2*(i-1)+1), '-', 'linewidth', linewith, 'Color', [0 0 0]/255); hold on;
%     xlim([0 max(T)]);
%     ylabel(xlb{i});
% end
% xlabel('sample time k');
% legend('Reference','H1', 'H10');

figure(5);
for i=1:2
    subplot(2,1,i);
    j=2*(i-1)+1;
    b1 =bar(T, reshape(obj2.Smat(j,j,T), 1, max(T)), 'hist');  hold on;
    set(b1, 'FaceColor', [0 0 255]/255, 'EdgeColor', [0 0 0]/255);
    b2 =bar(T, -reshape(obj3.Smat(j,j,T), 1, max(T)), 'hist'); hold on;
    set(b2, 'FaceColor', [255 0 0]/255, 'EdgeColor', [0 0 0]/255);
    xlim([0 max(T)]);
    ylabel(xlb{i});
    if i==1
        title({'Propagated Variances', titletex1});
    end
end
set(gcf, 'Renderer', 'painters');  % fix the shades problem of bar variances
xlabel('sample time k');
l1 = legend('GPMPC1', 'GPMPC2', 'Location', 'best');
set(l1,'FontName','Helvetica', 'FontSize', 10, 'LineWidth', 1);
% % myzoom;
% figure(4);
% for i=1:2
%     subplot(2,1,i);
%     plot(lgpobj1.u(i,:)', 'k-', 'linewidth', linewith);  hold on;
%     plot(lgpobj2.u(i,:)', 'b-.', 'linewidth', linewith); 
%     xlim([0 max(T)]);
%     ylabel(xlb{i});
% end
% xlabel('sample time k');
% legend('H1', 'H10');
% % myzoom;

% err1 = obj1.path(T, :)-obj1.x(T,[1,3]);
% err1 = sum(err1.^2)/length(T)
% 
% err2 = obj1.path(T, :)-obj2.gppred(T,[1,3]);
% err2 = sum(err2.^2)/length(T)
% 
% err3 = obj1.path(T, :)-obj3.gppred(T,[1,3]);
% err3 = sum(err3.^2)/length(T)
end

function obj = wrapObj(temp)
        obj.refpred = temp.refpred;
        obj.gppred = temp.gppred;
        obj.u = temp.input;
        obj.fval = temp.fval;
        obj.Smat = temp.Smat;
end

% % save eps with 'Arial' Font
% %----------------------------------------------------------
% printeps(1,strcat('pmpc_', pathtype, '_output')); 
% printeps(2,strcat('pmpc_', pathtype, '_predVar')); 
% printeps(3,strcat('pmpc_', pathtype, '_u1u2')); 