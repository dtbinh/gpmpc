clear; clc; close all;

%%  unfinished 
%%
pathtype = 'Lorenz';   % 'Lorenz'  'Helix'
switch pathtype
    case 'Elliptical';
        temp = load('Tobj_fminsearch_ell.mat'); Tobj1 = temp.obj;
        temp = load('Robj_fminsearch_ell.mat'); Robj1 = temp.obj;
        temp = load('Tobj_minimize_ell.mat'); Tobj2 = temp.obj;
        temp = load('Robj_minimize_ell.mat'); Robj2 = temp.obj;
        temp = load('Tobj_ltvgpmpc_ell.mat'); Tobj3 = temp.Tobj;
        temp = load('Robj_ltvgpmpc_ell.mat'); Robj3 = temp.Robj;
        T=1:189;
        path = pathGenerate(0.5, 200+1, pathtype);
        trajectory = [path.xr(T)',path.yr(T)',path.zr(T)'];
    case 'Lorenz';
        temp = load('Tobj_fminsearch_lorenz.mat');Tobj1 = temp.obj;
        temp = load('Robj_fminsearch_lorenz.mat'); Robj1 = temp.obj
        temp = load('Tobj_minimize_lorenz.mat');Tobj2 = temp.obj;
        temp = load('Robj_minimize_lorenz.mat'); Robj2 = temp.obj;
        temp = load('Tobj_ltvgpmpc_lorenz.mat'); Tobj3 = temp.Tobj;
        temp = load('Robj_ltvgpmpc_lorenz.mat'); Robj3 = temp.Robj;
        trajectory = Tobj2.refpred(1:170,:);
        T=1:170;
    case 'Helix';
        temp = load('Tobj_fminsearch_helix.mat');Tobj1 = temp.obj;
        temp = load('Robj_fminsearch_helix.mat'); Robj1 = temp.obj
        temp = load('Tobj_minimize_helix.mat');Tobj2 = temp.obj;
        temp = load('Robj_minimize_helix.mat'); Robj2 = temp.obj;
        temp = load('Tobj_ltvgpmpc_helix.mat'); Tobj3 = temp.Tobj;
        temp = load('Robj_ltvgpmpc_helix.mat'); Robj3 = temp.Robj;
        T=1:189;
        path = pathGenerate(0.5, 200+1, pathtype);
        trajectory = [path.xr(T)',path.yr(T)',path.zr(T)'];
end

linewith = 1.5;
figure(1);
xlb = {'X', 'Y', 'Z'};
for i=1:3
    subplot(3,1,i);
    plot(Tobj1.refpred(:,i),':', 'linewidth', linewith,'Color', [255 0 0]/255); hold on;
    plot(Tobj1.gppred(:,i), 'b-.', 'linewidth', 1);  hold on;
    plot(Tobj2.gppred(:,i), '-', 'linewidth', 1, 'Color', [105 105 105]/255); hold on;
    plot(Tobj3.gppred(:,i), '-', 'linewidth', 1, 'Color', [0 0 0]/255); 
    xlim([0 max(T)]);
    ylabel(xlb{i});
end
xlabel('Time(s)');
legend('Reference','PMPC', 'G-PMPC', 'LL-PMPC');

% figure(2);
% for i=1:3
%     subplot(3,1,i);
%     j=2*(i-1)+1;
%     b1 =bar(T, reshape(Tobj1.Smat(j,j,T), max(T),1), 1);  hold on;
%     set(b1,'FaceColor', [0 0 255]/255, 'EdgeColor', [0 0 0]/255);
%     set(get(b1,'Children'),'FaceAlpha',1);
%     b2 =bar(T, reshape(Tobj2.Smat(j,j,T), max(T),1), 1); hold on;
%     set(b2,'FaceColor', [255 0 0]/255, 'EdgeColor', [0 0 0]/255);
%     set(get(b2,'Children'),'FaceAlpha',1);
%     b3 =bar(T, reshape(Tobj3.S1(j,j,T), max(T),1), 1); hold on;
%     set(b3,'FaceColor', 'none', 'EdgeColor', [0 0 0]/255);
%     xlim([0 max(T)]);
%     ylabel(xlb{i});
% end
% xlabel('Times');
% l1 = legend('PMPC', 'G-PMPC', 'LL-PMPC', 'Location', 'best');
% set(l1,'FontName','Helvetica', 'FontSize', 10, 'LineWidth', 1);
% % myzoom;

figure(2);
for i=1:3
    j=2*(i-1)+1;
    subplot(3,1,i);
    barinfo = [reshape(Tobj1.Smat(j,j,T), max(T),1), ...
        reshape(Tobj2.Smat(j,j,T), max(T),1)];
    b1 =bar(T, -barinfo, 1);  hold on;
    set(b1(1),'FaceColor', [0 0 255]/255, 'EdgeColor', [0 0 0]/255);
    set(b1(2),'FaceColor', [255 0 0]/255, 'EdgeColor', [0 0 0]/255);
%     set(b1,'FaceColor', [0 0 255]/255, 'EdgeColor', [0 0 0]/255);
%     set(get(b1,'Children'),'FaceAlpha',1);
    b3 =bar(T, reshape(Tobj3.S1(j,j,T), max(T),1), 1); hold on;
    set(b3,'FaceColor', 'none', 'EdgeColor', [0 0 0]/255);
    xlim([0 max(T)]);
    ylabel(xlb{i});
end
xlabel('Time(s)');
l1 = legend('PMPC', 'G-PMPC', 'LL-PMPC', 'Location', 'best');
set(l1,'FontName','Helvetica', 'FontSize', 10, 'LineWidth', 1);
% myzoom;

figure(3);
xlb={'u_1', 'u_x', 'u_y'};
for i=1:3
    subplot(3,1,i);
    stairs(Tobj1.u(:,i), 'r--', 'linewidth', linewith); hold on;
    stairs(Tobj2.u(:,i), 'b-.', 'linewidth', linewith);  hold on;
    stairs(Tobj3.u(i,:)', 'k-', 'linewidth', 1); 
    xlim([0 max(T)]);
    ylabel(xlb{i});
end
xlabel('Time(s)');
legend('PMPC', 'G-PMPC','LL-PMPC');

figure(4);
xlb = {'\phi', '\theta', '\psi'};
for i=1:3
    subplot(3,1,i);
    plot(Robj1.refpred(:,i),':', 'linewidth', linewith,'Color', [255 0 0]/255); hold on;
    plot(Robj1.gppred(:,i), 'b-.', 'linewidth', 1);  hold on;
    plot(Robj2.gppred(:,i), '-', 'linewidth', 1, 'Color', [105 105 105]/255); hold on;
    plot(Robj3.gppred(:,i), '-', 'linewidth', 1, 'Color', [0 0 0]/255); 
    xlim([0 max(T)]);
    ylabel(xlb{i});
end
xlabel('Time(s)');
legend('Reference','PMPC', 'G-PMPC', 'LL-PMPC');

figure(5);
for i=1:3
    j=2*(i-1)+1;
    subplot(3,1,i);
    barinfo = [reshape(Robj1.Smat(j,j,T), max(T),1), ...
        reshape(Robj2.Smat(j,j,T), max(T),1)];
    b1 =bar(T, -barinfo, 1);  hold on;
    set(b1(1),'FaceColor', [0 0 255]/255, 'EdgeColor', [0 0 0]/255);
    set(b1(2),'FaceColor', [255 0 0]/255, 'EdgeColor', [0 0 0]/255);
%     set(b1,'FaceColor', [0 0 255]/255, 'EdgeColor', [0 0 0]/255);
%     set(get(b1,'Children'),'FaceAlpha',1);
    b3 =bar(T, reshape(Robj3.S1(j,j,T), max(T),1), 1); hold on;
    set(b3,'FaceColor', 'none', 'EdgeColor', [0 0 0]/255);
    xlim([0 max(T)]);
    ylabel(xlb{i});
end
xlabel('Time(s)');
l1 = legend('PMPC', 'G-PMPC', 'LL-PMPC', 'Location', 'best');
set(l1,'FontName','Helvetica', 'FontSize', 10, 'LineWidth', 1);
% myzoom;

figure(6);
xlb={'u_2', 'u_3', 'u_4'};
for i=1:3
    subplot(3,1,i);
    stairs(Robj1.u(:,i), 'r--', 'linewidth', linewith); hold on;
    stairs(Robj2.u(:,i), 'b-.', 'linewidth', linewith);  hold on;
    stairs(Robj3.u(i,:)', 'k-', 'linewidth', 1); 
    xlim([0 max(T)]);
    ylabel(xlb{i});
end
xlabel('Time(s)');
legend('PMPC', 'G-PMPC','LL-PMPC');

paraset = getRouteParams(pathtype);
quadTranRotaplot3(T, Tobj1.gppred, trajectory, ...
            Robj1.gppred, 'none', [], ...
            paraset.flyerpara,paraset.viewangle);

