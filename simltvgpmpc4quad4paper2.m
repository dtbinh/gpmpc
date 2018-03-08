function simltvgpmpc4quad4paper2

clear; clc; close all;
set(0, 'DefaultAxesFontName', 'Arial');  
set(0, 'DefaultTextFontName', 'Arial');  
pathtype = 'Elliptical';   % 'Lorenz'  'Elliptical'

switch pathtype
    case 'Elliptical';
        obj1 = load('gpmpc_quad_elliptical_quadprog.mat'); 
        obj2 = load('smpc_quad_Elliptical_MS.mat');
    case 'Lorenz';
        obj1 = load('gpmpc_quad_lorenz_quadprog.mat');
        obj2 = load('smpc_quad_Lorenz_MS.mat');
end

% quadTranRotaplot3(obj1.plotvec_t, obj1.plotvec_y, obj1.plotvec_ref, ...
%             obj1.plotvec_angle, 'none', [], ...
%             obj1.paraset.flyerpara,obj1.paraset.viewangle); hold on;

path2follow1 = strcat('ds_nmpc_', pathtype, '_trans', '.txt');
path2follow2 = strcat('ds_nmpc_', pathtype, '_rotate', '.txt');
temp1 = load(path2follow1);
temp2 = load(path2follow2);

figure(2);
linew = 1;
ylb = {'X[m]','Y[m]','Z[m]'};
for i=1:3
    subplot(3,1,i);
    plot(obj1.plotvec_ref(:,i), 'r:','linewidth', 1.5); hold on;
    plot(temp1(2:end, 2*(i-1)+1), '--','linewidth', linew, 'Color', [0.4 0.4 0.4]); hold on;
    plot(obj2.plotvec_y(:,i), 'k-.', 'linewidth', linew); hold on;
    plot(obj1.plotvec_y(:,i), 'b-', 'linewidth', linew);  
    ylabel(ylb{i});
    xlim([0 189]);
    if i==1
        title('Translational Subsystem Tracking');
    else
    end
end
legend('Reference', 'NMPC', 'Nonlinear GPMPC', 'GPMPC');
xlabel('sample time k');

iae1 = abs(obj1.plotvec_y-temp1(2:190,[1 3 5]));
iae2 = abs(obj2.plotvec_y-temp1(2:190,[1 3 5]));
figure(3);
for i=1:3
    subplot(3,1,i);
    plot(cumsum(iae2(:,i)), 'b--', 'linewidth', linew); hold on;
    plot(cumsum(iae1(:,i)), 'k-', 'linewidth', linew); hold on; 
    xlim([0 189]);
    ylabel(ylb{i});
    if i==1
        title('Integral Absolute Error (IAE) -- Translational Subsystem');
    else
    end
end
legend('Nonlinear GPMPC', 'GPMPC');
xlabel('sample time k');

plotvec_angleref = obj1.r_gppred;
figure(4);
ylb = {'\phi[rad]','\theta[rad]','\psi[rad]'};
for i=1:3
    subplot(3,1,i);
    plot(plotvec_angleref(:,i), 'r:','linewidth', 1.5);hold on;
    plot(temp2(2:end, 2*(i-1)+1), '--','linewidth', linew, 'Color', [0.4 0.4 0.4]); hold on;
    plot(obj2.plotvec_angle(:,i), 'k-.', 'linewidth', linew);
    plot(obj1.plotvec_angle(:,i), 'b-', 'linewidth', linew);
    ylabel(ylb{i});
    xlim([0 189]);
    if i==1
        title('Rotational Subsystem Tracking');
    else
    end
end
legend('Reference', 'NMPC', 'Nonlinear GPMPC','GPMPC');
xlabel('sample time k');

iae1 = abs(obj1.plotvec_angle-temp2(2:190,[1 3 5]));
iae2 = abs(obj2.plotvec_angle-temp2(2:190,[1 3 5]));
figure(5);
for i=1:3
    subplot(3,1,i);
    plot(cumsum(iae2(:,i)), 'b--', 'linewidth', linew); hold on;
    plot(cumsum(iae1(:,i)), 'k-', 'linewidth', linew); hold on;
    xlim([0 189]);
    ylabel(ylb{i});
    if i==1
        title('Integral Absolute Error (IAE) -- Rotational Subsystem');
    else
    end
end
legend('Nonlinear GPMPC', 'GPMPC');
xlabel('sample time k');

figure(6);
% uaxis = {[0 189 -45 0], [0 189 -2 2], [0 189 -2 2]};
uaxis = {[0 189 1 100], [0 189 -0.2 0.2], [0 189 -0.2 0.2]};
controlinputs1 = obj1.t_u';
controlinputs2 = obj2.t_u';
ylb = {'U_1','u_x','u_y'};
for i=1:3
    subplot(3,1,i);
    plot(temp1(:, 6+i), 'r:','linewidth', linew); hold on;
    plot(controlinputs2(:,i), 'k--', 'linewidth', linew);  hold on;
    plot(controlinputs1(:,i), 'b-', 'linewidth', linew); 
    ylabel(ylb{i});
    xlim([0 189]);
    if i==1
        title('Control Inputs in Translational Subsystem');
    else
    end
end
legend('NMPC', 'Nonlinear GPMPC', 'GPMPC');
xlabel('sample time k');

figure(7);
controlinputs1 = obj1.r_u';
controlinputs2 = obj2.r_u';
ylb = {'U_2','U_3','U_4'};
for i=1:3
    subplot(3,1,i);
    plot(temp2(:, 6+i), 'r:','linewidth', linew);  hold on;
    plot(controlinputs2(:,i), 'k--', 'linewidth', linew);  hold on;
    plot(controlinputs1(:,i), 'b-', 'linewidth', linew); 
    ylabel(ylb{i});
%     axis([0 189 0.1 0.9]);
    xlim([0 189]);
    if i==1
        title('Control Inputs in Rotational Subsystem');
    else
    end
end
legend('NMPC', 'Nonlinear GPMPC', 'GPMPC');
xlabel('sample time k');

end