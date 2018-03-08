
clc; clear; close all;

syscase = 'nltv';     %     ltv  nltv
loadType = 'original';
optimizer = 'minimize';    % 'sqp'  'ltvmpc'  'quadprog'
pathcase = 'lorenz';

% fminsearch  fminlbfgs  fmincon   minimize
[gppred, refpred, input, fval, Smat] = ...
    pgp_ltvmpc_nlsys(syscase, optimizer, loadType, pathcase);

path2follow = strcat('ds_nmpc_nltv_', pathcase,'.txt');
temp = load(path2follow);
nmpc_out = temp(2:end, [1 3]);
nmpc_in = temp(:, [5 6]);

switch pathcase
    case 'duffing';
        titletex1 = '"Duffing" Trajectory Tracking';
        titletex2 = 'Control inputs in "Duffing" Trajectory Tracking';
    case 'lorenz';
        titletex1 = '"Lorenz" Trajectory Tracking';
        titletex2 = 'Control inputs in "Lorenz" Trajectory Tracking';
    case 'step';
        titletex1 = '"Step" Trajectory Tracking';
        titletex2 = 'Control inputs in "Step" Trajectory Tracking';
    case 'sincos';
        titletex1 = '"Curve" Trajectory Tracking';
        titletex2 = 'Control inputs in "Curve" Trajectory Tracking';
end

figure(1);
switch syscase
    case 'ltv';
        seq = [1 2];
    case 'nltv';
        seq = [1 3];
end
yyl = {'Y_1(k)', 'Y_2(k)'};
for i=1:length(seq)
    subplot(length(seq),1,i);
    plot(refpred(:,seq(i)), 'r:'); hold on;
    plot(nmpc_out(:, i), 'k--'); hold on;
    plot(gppred(:,seq(i)), 'b-'); hold on;
    ylabel(yyl{i});
    xlim([0 189]);
    if i==1
        title(titletex1);
    else
    end
end
xlabel('sampling time k');
legend('Reference', 'NMPC', 'GPMPC');

% figure(2);
% uu = input';
% uul = {'u_1(k)', 'u_2(k)'};
% for i=1:size(uu,2)
%     subplot(size(uu,2),1,i);
%     plot(nmpc_in(:,i), 'k--'); hold on;
%     plot(uu(2:end,i), 'b-'); 
%     ylabel(uul{i});
%     xlim([0 189]);
%     if i==1
%         title(titletex2);
%     end
% end
% xlabel('sampling time k');
% legend('NMPC', 'GPMPC');

% figure(3);
% T=1:size(Smat,3);
% for j=1:2
%     subplot(2,1,j);
%     b1 =  bar(T, reshape(Smat(j,j,T), max(T),1)); hold on;
% end

% figure(4);
% plot(fval); 
% ylabel('Cost Values'); xlabel('sampling time k');
