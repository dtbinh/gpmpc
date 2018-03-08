clc; clear; close all;

% Lorenz   Elliptical   Helix   crown   spraypool   randpath
pathtype = 'Elliptical'; 
mpciterations = 199;   %  180    199   89  120

% fminsearch  fminlbfgs  fmincon   minimize
[t_gppred, t_refpred, t_lsyspred, t_u, t_fval, t_S] = ...
    pgp_ltvmpc_quad(pathtype, mpciterations, 'TRAN');

[r_gppred, r_refpred, r_lsyspred, r_u, r_fval, r_S] = ...
    pgp_ltvmpc_quad(pathtype, mpciterations, 'ROTATE');

paraset = getRouteParams(pathtype);
plotvec_t = 1:mpciterations-10;
plotvec_y = t_gppred;
plotvec_ref = t_refpred(plotvec_t,:);
plotvec_angle = r_gppred;

% Tobj.gppred = t_gppred;
% Tobj.refpred = t_refpred;
% Tobj.u = t_u;
% Tobj.S1 = t_S;
% 
% Robj.gppred = r_gppred;
% Robj.refpred = r_refpred;
% Robj.u = r_u;
% Robj.S1 = r_S;

quadTranRotaplot3(plotvec_t, plotvec_y, plotvec_ref, ...
            plotvec_angle, 'none', [], ...
            paraset.flyerpara,paraset.viewangle);
        
% dynamicQuadPlot3(plotvec_t, plotvec_y, plotvec_ref, t_u');

path2follow1 = strcat('ds_nmpc_', pathtype, '_trans', '.txt');
path2follow2 = strcat('ds_nmpc_', pathtype, '_rotate', '.txt');
temp1 = load(path2follow1);
temp2 = load(path2follow2);

figure(2);
ylb = {'X[m]','Y[m]','Z[m]'};
linew = 1;
for i=1:3
    subplot(3,1,i);
    plot(plotvec_ref(:,i), 'r:','linewidth', linew); hold on;
    plot(temp1(2:end, 2*(i-1)+1), 'k--','linewidth', linew); hold on;
    plot(plotvec_y(:,i), 'b-', 'linewidth', linew); 
    ylabel(ylb{i});
    xlim([0 189]);
    if i==1
        title('Translational Subsystem Tracking');
    else
    end
end
legend('Reference', 'NMPC', 'GPMPC');
xlabel('sample time k');

plotvec_angleref = r_gppred;
figure(3);
ylb = {'\phi[rad]','\theta[rad]','\psi[rad]'};
for i=1:3
    subplot(3,1,i);
    plot(plotvec_angleref(:,i), 'r:','linewidth', linew);hold on;
    plot(temp2(2:end, 2*(i-1)+1), 'k--','linewidth', linew); hold on;
    plot(plotvec_angle(:,i), 'b-', 'linewidth', linew);
    ylabel(ylb{i});
    xlim([0 189]);
    if i==1
        title('Rotational Subsystem Tracking');
    else
    end
end
legend('Reference', 'NMPC', 'GPMPC');
xlabel('sample time k');

figure(4);
% uaxis = {[0 189 -45 0], [0 189 -2 2], [0 189 -2 2]};
uaxis = {[0 189 1 100], [0 189 -0.2 0.2], [0 189 -0.2 0.2]};
controlinputs = t_u';
ylb = {'U_1','u_x','u_y'};
for i=1:3
    subplot(3,1,i);
    plot(temp1(:, 6+i), 'k-','linewidth', linew); hold on;
    plot(controlinputs(:,i), 'b-', 'linewidth', linew); 
    ylabel(ylb{i});
    xlim([0 189]);
    if i==1
        title('Control Inputs in Translational Subsystem');
    else
    end
end
legend('NMPC', 'GPMPC');
xlabel('sample time k');

figure(5);
controlinputs = r_u';
ylb = {'U_2','U_3','U_4'};
for i=1:3
    subplot(3,1,i);
    plot(temp2(:, 6+i), 'k-','linewidth', linew);  hold on;
    plot(controlinputs(:,i), 'b-', 'linewidth', linew); 
    ylabel(ylb{i});
%     axis([0 189 0.1 0.9]);
    xlim([0 189]);
    if i==1
        title('Control Inputs in Rotational Subsystem');
    else
    end
end
legend('NMPC', 'GPMPC');
xlabel('sample time k');


% figure(2);
% drawRealPrediction(plotvec_ref,plotvec_y,0);
% legend('Trajectory', 'GP Mean');
% 
% figure(3);
% plotvec_ref = r_refpred(plotvec_t,:);
% drawRealPrediction(plotvec_ref,plotvec_angle,0);
% legend('Trajectory', 'GP Mean');
% 
% figure(4);
% cu = t_u';
% subplot(2,1,1); plot(cu(:,1), 'k-');  legend('u_1');
% axis([0 189 -45 0]);
% subplot(2,1,2); plot(cu(:,2), 'b-'); hold on;
% plot(cu(:,3), 'k-'); legend('u_x', 'u_y');
% axis([0 189 -2 2]);
% 
% figure(5);
% cu = r_u';
% subplot(2,1,1); plot(cu(:,1), 'k-'); hold on;
% plot(cu(:,2), 'b-'); legend('u_2', 'u_3');
% subplot(2,1,2); plot(cu(:,3), 'k-');legend( 'u_4');

% figure(6);
% subplot(1,2,1); plot(t_fval); 
% ylabel('Cost Values'); xlabel('Steps');
% subplot(1,2,2); plot(r_fval); xlabel('Steps');

% figure(7);
% T=1:mpciterations-10;
% for j=1:3
%     subplot(3,1,j);
%     b1 =  bar(T, reshape(t_S(j,j,T), max(T),1)); hold on;
% end


%  plot(t_gppred(:,1)); hold on; 
%  plot(t_gppred(:,1)+2*reshape(t_S(1,1,:), 189,1), 'k:'); hold on;
%  plot(t_gppred(:,1)-2*reshape(t_S(1,1,:), 189,1), 'k:'); hold on;
%  myzoom 
%  
% figure(7);
% i=2;
% xlb = {'X(m)', 'Y(m)', 'Z(m)'};
% subplot(2,1,1);
% plot(t_gppred(:,i), 'b-', 'markersize',2); hold on;
% plot(t_refpred(:,i),'k:','linewidth', 1);
% xlim([0 189]);
% ylabel(xlb{i});
% legend('Mean', 'Reference');
% subplot(2,1,2);
% bar(1:189,reshape(t_S(i,i,:), 189,1), 0.5);
% xlim([0 189]);
% xlabel('Times');
% ylabel('Confidence Intervals');
% myzoom;


