clear; clc; close all;
set(0, 'DefaultAxesFontName', 'Arial');  
set(0, 'DefaultTextFontName', 'Arial');  

pathtype = 'Lorenz';   % Lorenz    Elliptical
switch pathtype
    case 'Elliptical';
        temp = load('Tobj_fminsearch_ell.mat'); Tobj1 = temp.obj;
        temp = load('Robj_fminsearch_ell.mat'); Robj1 = temp.obj;
        temp = load('Tobj_minimize_ell.mat'); Tobj2 = temp.obj;
        temp = load('Robj_minimize_ell.mat'); Robj2 = temp.obj;
        T=1:189;
        path = pathGenerate(0.5, 200+1, pathtype);
        trajectory = [path.xr(T)',path.yr(T)',path.zr(T)'];
        temptxt='_ell';
    case 'Lorenz';
        temp = load('Tobj_fminsearch_lorenz.mat');Tobj1 = temp.obj;
        temp = load('Robj_fminsearch_lorenz.mat'); Robj1 = temp.obj;
        temp = load('Tobj_minimize_lorenz.mat');Tobj2 = temp.obj;
        temp = load('Robj_minimize_lorenz.mat'); Robj2 = temp.obj;
        trajectory = Tobj2.refpred(1:170,:);
        T=1:170;
        temptxt='_lorenz';
end

i=1;
switch i
    case {1,2};
        coe = 1;
    case 3;
        coe = 0;
end
temperr = Tobj1.refpred(T, i)-Tobj1.gppred(T,i);
gpmpc_tran_err = sum(temperr.^2)/length(T)
temperr = Tobj2.refpred(T, i)-Tobj2.gppred(T,i);
pmpc_tran_err = sum(temperr.^2)/length(T)

temperr = coe*Robj1.refpred(T, i)-Robj1.gppred(T,i);
gpmpc_rotate_err = sum(temperr.^2)/length(T)
temperr = coe*Robj2.refpred(T, i)-Robj2.gppred(T,i);
pmpc_rotate_err = sum(temperr.^2)/length(T)

linewith = 1.5;
figure(1);
xlb = {'X[m]', 'Y[m]', 'Z[m]'};
for i=1:3
    subplot(3,1,i);
    plot(Tobj1.refpred(:,i),'.-', 'linewidth', 1,'Color', [105 105 105]/255); hold on;
    plot(Tobj1.gppred(:,i), 'b-.', 'linewidth', linewith);  hold on;
    plot(Tobj2.gppred(:,i), 'k-', 'linewidth', 1); 
    xlim([0 max(T)]);
    ylabel(xlb{i});
end
xlabel('sample time k');
legend('Reference','GPMPC', 'Grad-GPMPC');

figure(2);
xlb = {'\phi [rad]', '\theta [rad]', '\psi [rad]'};
for i=1:2
    subplot(3,1,i);
    plot(Robj1.refpred(:,i),'.-', 'linewidth', 1,'Color', [105 105 105]/255); hold on;
    plot(Robj1.gppred(:,i), 'b-.', 'linewidth', linewith);  hold on;
    plot(Robj2.gppred(:,i), 'k-', 'linewidth', 1); 
    xlim([0 max(T)]);
    ylabel(xlb{i});
end
subplot(3,1,3);
plot(0.*Robj1.refpred(:,3),'.-', 'linewidth', 1,'Color', [105 105 105]/255); hold on;
plot(Robj1.gppred(:,3), 'b-.', 'linewidth', linewith);  hold on;
plot(Robj2.gppred(:,3), 'k-', 'linewidth', 1); 
xlim([0 max(T)]);
ylabel(xlb{3});
xlabel('sample time k');
legend('Reference', 'GPMPC', 'Grad-GPMPC');

% linewith = 1.5;
% figure(1);
% xlb = {'X[m]', 'Y[m]', 'Z[m]'};
% for i=1:3
%     subplot(3,1,i);
%     plot(Tobj1.refpred(:,i),':', 'linewidth', linewith,'Color', [255 0 0]/255); hold on;
%     plot(Tobj1.gppred(:,i), 'b-.', 'linewidth', linewith);  hold on;
%     plot(Tobj2.gppred(:,i), '-', 'linewidth', linewith, 'Color', [0 0 0]/255); 
%     xlim([0 max(T)]);
%     ylabel(xlb{i});
% end
% xlabel('sample time k');
% legend('Reference','PMPC', 'G-PMPC');
% 
% figure(2);
% xlb = {'variances on X', 'variances on Y', 'variances on Z'};
% for i=1:3
%     subplot(3,1,i);
%     j=2*(i-1)+1;
%     b1 =bar(T, reshape(Tobj1.Smat(j,j,T), max(T),1), 1);  hold on;
%     set(b1,'FaceColor', [0 0 255]/255, 'EdgeColor', [0 0 0]/255);
%     b2 =bar(T, -reshape(Tobj2.Smat(j,j,T), max(T),1), 1); hold on;
%     set(b2,'FaceColor', [255 0 0]/255, 'EdgeColor', [0 0 0]/255);
%     xlim([0 max(T)]);
%     ylabel(xlb{i});
% end
% xlabel('sample time k');
% set(gcf, 'Renderer', 'painters'); 
% l1 = legend('G-PGPMPC', 'PMPC', 'Location', 'best');
% set(l1,'FontName','Helvetica', 'FontSize', 10, 'LineWidth', 1);
% % myzoom;
% 
% figure(3);
% xlb={'u_1(k)', 'u_x(k)', 'u_y(k)'};
% for i=1:3
%     subplot(3,1,i);
%     stairs(Tobj1.u(:,i), 'k-', 'linewidth', linewith); hold on;
%     stairs(Tobj2.u(:,i), 'b-.', 'linewidth', linewith); 
%     xlim([0 max(T)]);
%     ylabel(xlb{i});
% end
% xlabel('sample time k');
% legend('PMPC', 'G-PMPC');
% 
% figure(4);
% xlb = {'\phi [rad]', '\theta [rad]', '\psi [rad]'};
% for i=1:2
%     subplot(3,1,i);
%     plot(Robj1.refpred(:,i),':', 'linewidth', linewith,'Color', [255 0 0]/255); hold on;
%     plot(Robj1.gppred(:,i), 'b-.', 'linewidth', linewith);  hold on;
%     plot(Robj2.gppred(:,i), '-', 'linewidth', linewith,'Color', [0 0 0]/255); 
%     xlim([0 max(T)]);
%     ylabel(xlb{i});
% end
% subplot(3,1,3);
% plot(0.*Robj1.refpred(:,3),':', 'linewidth', linewith,'Color', [255 0 0]/255); hold on;
% plot(Robj1.gppred(:,3), 'b-.', 'linewidth', linewith);  hold on;
% plot(Robj2.gppred(:,3), '-', 'linewidth', linewith,'Color', [0 0 0]/255); 
% xlim([0 max(T)]);
% ylabel(xlb{3});
% xlabel('sample time k');
% legend('Reference', 'PMPC', 'G-PMPC');
% 
% figure(5);
% xlb = {'variances on \phi', 'variances on \theta', 'variances on \psi'};
% for i=1:3
%     subplot(3,1,i);
%     j=2*(i-1)+1;
%     b1 =bar(T, reshape(Robj1.Smat(j,j,T), max(T),1), 1);  hold on;
%     set(b1,'FaceColor', [0 0 255]/255, 'EdgeColor', [0 0 0]/255);
%     b2 =bar(T, -reshape(Robj2.Smat(j,j,T), max(T),1), 1); hold on;
%     set(b2,'FaceColor', [255 0 0]/255, 'EdgeColor', [0 0 0]/255);
%     xlim([0 max(T)]);
%     ylabel(xlb{i});
% end
% xlabel('sample time k');
% l1 = legend('G-PMPC', 'PMPC', 'Location', 'best');
% set(l1,'FontName','Helvetica', 'FontSize', 10, 'LineWidth', 1);
% % myzoom;
% 
% figure(6);
% xlb={'u_2(k)', 'u_3(k)', 'u_4(k)'};
% for i=1:3
%     subplot(3,1,i);
%     stairs(Robj1.u(:,i), 'k-', 'linewidth', linewith); hold on;
%     stairs(Robj2.u(:,i), 'b-.', 'linewidth', linewith); 
%     xlim([0 max(T)]);
%     ylabel(xlb{i});
% end
% xlabel('sample time k');
% legend('PMPC', 'G-PMPC');
% 
% paraset = getRouteParams(pathtype);
% quadTranRotaplot3(T, Tobj1.gppred, trajectory, ...
%             Robj1.gppred, 'none', [], ...
%             paraset.flyerpara,paraset.viewangle);

% save eps with 'Arial' Font
%----------------------------------------------------------
% printeps(1,strcat('pmpc_', 'xyz', temptxt)); 
% printeps(2,strcat('pmpc_', 'xyzVariance', temptxt)); 
% printeps(3,strcat('pmpc_', 'xyzU', temptxt)); 
% printeps(4,strcat('pmpc_', 'angle', temptxt)); 
% printeps(5,strcat('pmpc_', 'angleVariance', temptxt)); 
% printeps(6,strcat('pmpc_', 'angleU', temptxt)); 
% printeps(7,strcat('pmpc_', 'track', temptxt)); 
