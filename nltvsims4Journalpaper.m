function nltvsims4Journalpaper
clear;clc;close all;
set(0, 'DefaultAxesFontName', 'Arial');  
set(0, 'DefaultTextFontName', 'Arial');  

pathtype = 'lorenz';     % duffing  lorenz  step
T=1:189;
switch pathtype
    case 'duffing';
        temp = load('nmpc_nltv_duffing.mat'); obj1 = temp.obj;
        temp = load('pmpc_duffing_sqp.mat'); obj2 = wrapObj(temp);
        temp = load('pmpc_duffing_lgpmpc.mat'); obj3 = wrapObj(temp);
    case 'lorenz';
        temp = load('nmpc_nltv_lorenz.mat'); obj1 = temp.obj;
        temp = load('pmpc_lorenz_sqp_H1.mat'); obj2 = wrapObj(temp);
        temp = load('pmpc_lorenz_sqp_H10.mat'); obj2= wrapObj(temp);
        temp = load('pmpc_lorenz_lgpmpc_H1.mat'); obj3 = wrapObj(temp);
        temp = load('pmpc_lorenz_lgpmpc_H10.mat'); obj3 = wrapObj(temp);
    case 'step';
        obj1.path =  TracjectorGenerator(1:200,4);
        temp = load('pmpc_step_sqp_H3.mat'); obj2 = wrapObj(temp);
        temp = load('pmpc_step_lgpmpc_H3.mat'); obj3 = wrapObj(temp);
end

linewith = 1.5;
xlb = {'Y_1(k)', 'Y_2(k)'};
figure(1);
for i=1:2
    subplot(2,1,i);
    plot(obj1.path(:,i),'k-'); hold on;
%     plot(obj1.x(:,2*(i-1)+1), 'g-', 'linewidth', linewith);  hold on;
    plot(obj2.gppred(:,2*(i-1)+1), 'b-', 'linewidth', linewith); hold on;
    plot(obj3.gppred(:,2*(i-1)+1), 'r:', 'linewidth', linewith); hold on;
%     plot(sqpobj3.gppred(:,2*(i-1)+1), '-', 'linewidth', linewith, 'Color', [0 0 0]/255); 
    xlim([0 max(T)]);
    ylabel(xlb{i});
end
xlabel('sample time k');
legend('Reference','GPMPC1', 'GPMPC2');

figure(2);
xlb={'u_1(k)', 'u_2(k)'};
for i=1:2
    subplot(2,1,i);
    plot(obj2.u(i,:), 'k--', 'linewidth', linewith); hold on;
    plot(obj3.u(i,:)', 'b-');  hold on;
%     plot(obj1.u(i,:)', 'r--', 'linewidth', linewith); 
    xlim([0 max(T)]);
    ylabel(xlb{i});
end
xlabel('sample time k');
legend('GPMPC1', 'GPMPC2');

iae1 = abs(obj2.gppred(:,[1 3])-obj1.path(1:189,:));
iae2 = abs(obj3.gppred(:,[1 3])-obj1.path(1:189,:));
figure(3);
subplot(2,1,1);
plot(cumsum(iae1(:,1)), 'k--', 'linewidth', linewith); hold on;
plot(cumsum(iae2(:,1)), 'b-'); 
xlim([0 max(T)]);
ylabel('IAE on Y_1');
subplot(2,1,2);
plot(cumsum(iae1(:,2)), 'k--', 'linewidth', linewith); hold on;
plot(cumsum(iae2(:,2)), 'b-'); 
xlim([0 max(T)]);
ylabel('IAE on Y_2');
xlabel('sample time k');
legend('GPMPC1', 'GPMPC2');

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

figure(4);
xlb = {'variances on Y_1(k)', 'variances on Y_2(k)'};
for i=1:2
    subplot(2,1,i);
    j=2*(i-1)+1;
    b1 =bar(T, reshape(obj2.Smat(j,j,T), 1, max(T)), 'hist');  hold on;
    set(b1, 'FaceColor', [0 0 255]/255, 'EdgeColor', [0 0 0]/255);
    b2 =bar(T, -reshape(obj3.Smat(j,j,T), 1, max(T)), 'hist'); hold on;
    set(b2, 'FaceColor', [255 0 0]/255, 'EdgeColor', [0 0 0]/255);
    xlim([0 max(T)]);
    ylabel(xlb{i});
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




