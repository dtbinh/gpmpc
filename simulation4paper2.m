clear;clc;close all;
ttype = 'ltv2';

switch ttype
    case 'ltv1';
        obj1= load('ltvsysobj.mat');
        obj2= load('ltvobj_fminsearch.mat');
        obj3= load('ltvobj_minimize.mat');
        T = 1:250;
    case 'ltv2';
        obj1= load('lmpc4ltv_lorenz.mat');
        obj2= load('ltvobj_lorenz_fminsearch.mat');
        obj3= load('ltvobj_lorenz_minimize.mat');
        T = 1:180;
end

obj1 = obj1.obj;
obj2 = obj2.obj;
obj3 = obj3.obj;

i=1;
subplot(2,1,1);
plot(obj1.path(T,i), '.-', 'linewidth', 1.5,'Color', [105 105 105]/255); hold on;
plot(obj2.refpred(:,i),':', 'linewidth', 1.5,'Color', [255 0 0]/255); hold on;
plot(obj2.gppred(:,i), 'b-.', 'linewidth', 1.5); hold on;
plot(obj3.gppred(:,i), '-', 'linewidth', 1.5,'Color', [0 0 0]/255); 
ylabel('x_1(k)');xlim([0 max(T)]);
i=2;
subplot(2,1,2);
plot(obj1.path(T,i), '.-', 'linewidth', 1.5,'Color', [105 105 105]/255); hold on;
plot(obj2.refpred(:,i),':', 'linewidth', 1.5,'Color', [255 0 0]/255); hold on;
plot(obj2.gppred(:,i), 'b-.', 'linewidth', 1.5); hold on;
plot(obj3.gppred(:,i), '-', 'linewidth', 1.5,'Color', [0 0 0]/255); 
ylabel('x_2(k)');xlim([0 max(T)]);
xlabel('sample time k');
legend('Reference','LMPC', 'GPMPC', 'Grad-based GPMPC');

i=1;
figure(2);
subplot(2,1,1);
b1 =bar(T, reshape(obj3.Smat(i,i,:), max(T),1), 1);  hold on;
set(b1,'FaceColor', [255 0 0]/255, 'EdgeColor', [0 0 0]/255);
b2 =bar(T, -reshape(obj2.Smat(i,i,:), max(T),1), 1);
set(b2,'FaceColor', [0 0 255]/255, 'EdgeColor', [0 0 0]/255);
xlim([0 max(T)]);ylabel('variances on x_1');
i=2;
subplot(2,1,2);
b1 =bar(T, reshape(obj3.Smat(i,i,:), max(T),1), 1);  hold on;
set(b1,'FaceColor', [255 0 0]/255, 'EdgeColor', [0 0 0]/255);
b2 =bar(T, -reshape(obj2.Smat(i,i,:), max(T),1), 1);
set(b2,'FaceColor', [0 0 255]/255, 'EdgeColor', [0 0 0]/255);
xlim([0 max(T)]);
ylabel('variances on x_2');xlabel('sample time k');
legend('Grad-based GPMPC','GPMPC');
% myzoom;

i=1;
figure(3);subplot(2,1,1);
plot(obj1.u(:,i),':', 'linewidth', 1.5,'Color', [255 0 0]/255); hold on;
plot(obj2.u(i,:)', 'b-.', 'linewidth', 1.5); hold on;
plot(obj3.u(i,:)', '-', 'linewidth', 1.5,'Color', [0 0 0]/255);
xlim([0 max(T)]);ylabel('u_1(k)');
i=2;subplot(2,1,2);
plot(obj1.u(:,i),':', 'linewidth', 1.5,'Color', [255 0 0]/255); hold on;
plot(obj2.u(i,:)', 'b-.', 'linewidth', 1.5); hold on;
plot(obj3.u(i,:)', '-', 'linewidth', 1.5,'Color', [0 0 0]/255);
legend('LMPC', 'GPMPC', 'Grad-based GPMPC');
xlim([0 max(T)]);
xlabel('sample time k');
ylabel('u_2(k)');

% figure(4);
% plot(obj1.path(:,1), obj1.path(:,2), 'k:', 'linewidth', 1.5); hold on;
% plot(obj2.gppred(:,1), obj2.gppred(:,2), 'r--', 'linewidth', 1.5); hold on;
% plot(obj3.gppred(:,1), obj3.gppred(:,2), 'b-', 'linewidth', 1.5); hold on;
% xlabel('Y_1');ylabel('Y_2');
% legend('Reference', 'GPMPC', 'Grad-based GPMPC');
% plot(obj1.path(1,1), obj1.path(1,2), 'k.', 'MarkerSize', 10,'linewidth', 1.5); hold on;
% plot(obj2.gppred(1,1), obj2.gppred(1,2), 'ro', 'MarkerSize', 10,'linewidth', 1.5); hold on;
% plot(obj2.gppred(1,1), obj2.gppred(1,2), 'b^', 'MarkerSize', 10,'linewidth', 1.5); hold on;

% computing tracking errors
% T = 25:250;
err1 = obj1.path(T, :)-obj1.x(T,:);
err1 = sum(err1.^2)/length(T);

err2 = obj1.path(T, :)-obj2.gppred(T,:);
err2 = sum(err2.^2)/length(T);

err3 = obj1.path(T, :)-obj3.gppred(T,:);
err3 = sum(err3.^2)/length(T);

