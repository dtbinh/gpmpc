%  NN based unknown NLTV system learning
clear; clc; close all;

nnSolver = 'rbf';   %   rbf  nn
normalornot = 'no';
samples = load('ds_nmpc_nltv_duffing.txt');
%   ds_nmpc_nltv_step
%   ds_nmpc_nltv_lorenz
%   ds_nmpc_nltv_sincos
%   ds_nmpc_nltv_duffing

dsx = samples(:,1:6);
dsy = samples(:,7:10);
obnum = 189;

switch normalornot
    case 'yes';
        [ddsx, Psx]=mapminmax(dsx'); 
        [ddsy, Psy]=mapminmax(dsy');
        xtrain=ddsx(:,1:obnum);
        ytrain=ddsy(:,1:obnum);
    case 'no';
        xtrain = dsx';
        ytrain = dsy';
        ddsx = dsx';
end


tic
switch nnSolver
    case 'nn';
        net = newff(xtrain, ytrain, [20],{'tansig'},'traincgb', 'learngd');
        net.trainParam.epochs=1000; 
        net.trainParam.goal= 1e-06;
        net.trainParam.showWindow= 0;
        net=train(net,xtrain,ytrain);
    case 'rbf';
        net = newrb(xtrain, ytrain, 0.0012);
end
toc

xtest = ddsx;
outputs = net(xtest);

switch normalornot
    case 'yes';
        revpre = mapminmax('reverse',outputs,Psy)';
        revtest = mapminmax('reverse',ddsy,Psy)';
        revpre = revpre(:,[1 3]);
        revtest = revtest(:,[1 3]);
    case 'no';
        outputs = outputs';
        revpre = outputs(:,[1 3]);
        revtest = ytrain([1 3], :)';
end

figure(1);
subplot(2,1,1);
plot(revtest(:,1), 'r:'); hold on;
plot(revpre(:,1), 'b-'); 
ylabel('Y_1(k)');
subplot(2,1,2);
plot(revtest(:,2), 'r:'); hold on;
plot(revpre(:,2), 'b-'); 
ylabel('Y_2(k)');

figure(2);
subplot(2,1,1);
plot(revtest(:,1)-revpre(:,1), 'k', 'LineWidth', 1.5);
ylabel('Y_1(k)');
subplot(2,1,2);
plot(revtest(:,2)-revpre(:,2), 'k', 'LineWidth', 1.5);
ylabel('Y_2(k)');
xlabel('sample time k');
xlim([0 189]);

trainMSE = sum(sum((revpre-revtest).^2))/obnum


 
