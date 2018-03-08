clear;clc;close all;
set(0, 'DefaultAxesFontName', 'Arial');  
set(0, 'DefaultTextFontName', 'Arial');  

pathtype = 'duffing';   %  sincos   duffing  lorenz
switch pathtype
    case 'sincos';
        temp = load('nmpc_nltv_sincos.mat'); obj1 = temp.obj;
        temp = load('pmpc_sincos_fminsearch1.mat'); obj2 = temp.obj;
        temp = load('pmpc_sincos_minimize1.mat'); obj3 = temp.obj;
        T=1:160;
    case 'duffing';
        temp = load('nmpc_nltv_duffing.mat'); obj1 = temp.obj;
        temp = load('pmpc_duffing_fminsearch.mat'); obj2 = temp.obj;
        temp = load('pmpc_duffing_minimize.mat'); obj3 = temp.obj;
        T=1:189;
    case 'lorenz';
        temp = load('nmpc_nltv_lorenz.mat'); obj1 = temp.obj;
        temp = load('pmpc_lorenz_fminsearch.mat'); obj2 = temp.obj;
        temp = load('pmpc_lorenz_minimize.mat'); obj3 = temp.obj;
        T=1:189;
end

linewith = 1.5;
figure(1);
xlb = {'y_1(k)', 'y_2(k)'};
for i=1:2
    subplot(2,1,i);
    plot(obj1.path(:,i),':', 'linewidth', linewith,'Color', [255 0 0]/255); hold on;
    plot(obj1.x(:,2*(i-1)+1), 'b-.', 'linewidth', linewith);  hold on;
    plot(obj2.gppred(:,2*(i-1)+1), '-', 'linewidth', linewith, 'Color', [105 105 105]/255); hold on;
    plot(obj3.gppred(:,2*(i-1)+1), '-', 'linewidth', linewith, 'Color', [0 0 0]/255); 
    xlim([0 max(T)]);
    ylabel(xlb{i});
end
xlabel('sample time k');
legend('Reference','NMPC', 'GPMPC', 'Grad-based GPMPC');

figure(2);
xlb = {'variances on y_1(k)', 'variances on y_2(k)'};
linewith =1;
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
l1 = legend('Grad-based GPMPC', 'GPMPC', 'Location', 'best');
set(l1,'FontName','Helvetica', 'FontSize', 10, 'LineWidth', 1);
% myzoom;

linewith =1.5;
figure(3);
xlb={'u_1(k)', 'u_2(k)'};
for i=1:2
    subplot(2,1,i);
    stairs(obj1.u(:,i), 'r--', 'linewidth', linewith); hold on;
    stairs(obj2.u(i,:)', 'b-.', 'linewidth', linewith);  hold on;
    stairs(obj3.u(i,:)', 'k-', 'linewidth', linewith); 
    xlim([0 max(T)]);
    ylabel(xlb{i});
end
xlabel('sample time k');
legend('NMPC', 'GPMPC','Grad-based GPMPC');

err1 = obj1.path(T, :)-obj1.x(T,[1,3]);
err1 = sum(err1.^2)/length(T);

err2 = obj1.path(T, :)-obj2.gppred(T,[1,3]);
err2 = sum(err2.^2)/length(T);

err3 = obj1.path(T, :)-obj3.gppred(T,[1,3]);
err3 = sum(err3.^2)/length(T);

% % save eps with 'Arial' Font
% %----------------------------------------------------------
% printeps(1,strcat('pmpc_', pathtype, '_output')); 
% printeps(2,strcat('pmpc_', pathtype, '_predVar')); 
% printeps(3,strcat('pmpc_', pathtype, '_u1u2')); 




