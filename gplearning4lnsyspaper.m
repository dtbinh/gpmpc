clear; clc; close all;
set(0, 'DefaultAxesFontName', 'Arial');  
set(0, 'DefaultTextFontName', 'Arial');  

systype = 'nltv';   %  ltv  nltv
pathtype = 'lorenz';

switch systype
    case 'ltv';
        opt.state = 1:4;
        opt.in = 3:4;
        opt.out = 5:6;
        dynmodel.iNum =2;
        dynmodel.sNum =2;
        dynmodel.oNum =2; 
        
        switch pathtype
            case 'lorenz';
                filename = 'lmpc4ltv_lorenz';
                mpciterations = 180; %  180    199   170  250
                trainsize = 100;
                titletxt1 = {'"Lorenz" Trajectory -- Traning Errors', 'Traning Size = 100'};
                titletxt2 = {'"Lorenz" Trajectory -- Testing Errors', 'Testing Size = 180'};
            case 'duffing';
                filename = 'lmpc4ltvsys';
                mpciterations = 250; 
                trainsize = 140;
                titletxt1 = {'"Duffing" Trajectory -- Traning Errors', 'Traning Size = 140'};
                titletxt2 = {'"Duffing" Trajectory -- Testing Errors', 'Testing Size = 250'};
        end
    case 'nltv';
        mpciterations = 189;
        trainsize = 189;
        opt.state = 1:6;
        opt.in = 5:6;
        opt.out = 7:10;
%         filename = 'ds_nmpc_nltv_sincos';
%         filename = 'ds_nmpc_nltv_lorenz';
        filename = 'ds_nmpc_nltv_step';
%         filename = 'ds_nmpc_nltv_duffing';
        dynmodel.iNum =2;
        dynmodel.sNum =4;
        dynmodel.oNum =4; 
        titletxt1 = {'"Step" Trajectory -- Traning Errors', 'Traning Size = 189'};
        titletxt2 = {'"Step" Trajectory -- Testing Errors', 'Testing Size = 189'};
end

loadType = 'original';
opt.dsSampleRange = 1:trainsize;
opt.dsTestRange = 1:mpciterations;
trainds = loadSampleData(filename,loadType,opt);

dynmodel.inputs = trainds.proData(opt.dsSampleRange,opt.state);
dynmodel.targets = trainds.proData(opt.dsSampleRange,opt.out);
dynmodel.train = @train;


tic
plant = [];
[dynmodel nlml] = dynmodel.train(dynmodel, plant, 1000);
toc   % counting dynamic model learning time

traindata.inputs = trainds.proData(opt.dsSampleRange,opt.state);
traindata.targets = trainds.proData(opt.dsSampleRange,opt.out);

testdata.inputs = trainds.proData(opt.dsTestRange,opt.state);
testdata.targets = trainds.proData(opt.dsTestRange,opt.out);

S0size = dynmodel.sNum + dynmodel.iNum;
Ssize = dynmodel.sNum;

s= diag([0.1*ones(1, Ssize)].^2);
s0=  zeros(S0size,S0size,trainsize);
for i=1:trainsize
    s0(1:Ssize,1:Ssize,i)=s;
    [trainpred(i,:) s]= gp0d(dynmodel, traindata.inputs(i,:)', s0(:,:,i));
end
trainMSE = sum(sum((trainpred-traindata.targets).^2))/trainsize

figure(1);
subplot(2,1,1);
plot(trainpred(:,1)-traindata.targets(:,1), 'k', 'LineWidth', 1.5);
ylabel('Y_1(k)');
title(titletxt1);
xlim([0 trainsize]);
% xlim([0 190]);
% axis([0 190 -0.1 0.1]);
subplot(2,1,2);
plot(trainpred(:,2)-traindata.targets(:,2), 'k', 'LineWidth', 1.5);
ylabel('Y_2(k)');
xlabel('sample time k');
xlim([0 trainsize]);
% xlim([0 190]);
% axis([0 190 -0.1 0.1]);

% figure(1);
% subplot(4,1,1); i=1; 
% bar(1:trainsize, reshape(s0(i,i,:), trainsize,1), 1); hold on;
% xlim([0 trainsize]);
% subplot(4,1,2); i=2; 
% bar(1:trainsize, reshape(s0(i,i,:), trainsize,1), 1); hold on;
% xlim([0 trainsize]);
% subplot(4,1,3); i=3; 
% bar(1:trainsize, reshape(s0(i,i,:), trainsize,1), 1); hold on;
% xlim([0 trainsize]);
% subplot(4,1,4); i=4;
% bar(1:trainsize, reshape(s0(i,i,:), trainsize,1), 1); 
% xlim([0 trainsize]);

% s= diag([0.1*ones(1, Ssize)].^2);
% s1=  zeros(S0size,S0size,mpciterations);
% for i=1:mpciterations
%     s1(1:Ssize,1:Ssize,i)=s;
%     [testpred(i,:) s]= gp0d(dynmodel, testdata.inputs(i,:)', s1(:,:,i));
% end
% testMSE = sum(sum((testpred-testdata.targets).^2))/mpciterations
% 
% figure(2);
% subplot(2,1,1);
% plot(testpred(:,1)-testdata.targets(:,1), 'k', 'LineWidth', 1.5);
% ylabel('Y_1(k)');
% xlim([0 mpciterations]);
% title(titletxt2);
% axis([0 190 -0.1 0.1]);
% subplot(2,1,2);
% plot(testpred(:,2)-testdata.targets(:,2), 'k', 'LineWidth', 1.5);
% ylabel('Y_2(k)');
% xlim([0 mpciterations]);
% xlabel('sample time k');

% figure(1);
% i=1;
% b1=bar(1:mpciterations, reshape(s1(i,i,:), mpciterations,1), 1); hold on;
% set(b1,'FaceColor', [0 0 255]/255, 'EdgeColor', [0 0 0]/255);
% axis([0 mpciterations 0 0.02]);
% xlabel('Times');
% ylabel('x_1');

% figure(2);
% subplot(4,1,1); 
% i=1;
% b1=bar(1:mpciterations, reshape(s1(i,i,:), mpciterations,1), 1); hold on;
% set(b1,'FaceColor', [0 0 255]/255, 'EdgeColor', [0 0 0]/255);
% % axis([0 mpciterations 0 0.015]);
% ylabel('x_1');
% subplot(4,1,2); i=2; 
% b2=bar(1:mpciterations, reshape(s1(i,i,:), mpciterations,1), 1); hold on;
% set(b2,'FaceColor', [0 0 255]/255, 'EdgeColor', [0 0 0]/255);
% ylabel('x_2');
% axis([0 mpciterations 0 0.015]);
% subplot(4,1,3); i=3; 
% b3=bar(1:mpciterations, reshape(s1(i,i,:), mpciterations,1), 1); hold on;
% set(b3,'FaceColor', [0 0 255]/255, 'EdgeColor', [0 0 0]/255);
% ylabel('x_3');
% axis([0 mpciterations 0 0.015]);
% subplot(4,1,4); i=4;
% b4=bar(1:mpciterations, reshape(s1(i,i,:), mpciterations,1), 1); 
% set(b4,'FaceColor', [0 0 255]/255, 'EdgeColor', [0 0 0]/255);
% axis([0 mpciterations 0 0.015]);
% ylabel('x_4');
% xlabel('Time(s)');

% figure(2);
% xlb = {'variances on x_1', 'variances on x_2'};
% for i=1:2
%     subplot(2,1,i);
%     b=bar(1:mpciterations, ...
%         reshape(s1(2*(i-1)+1,2*(i-1)+1,:), ...
%         mpciterations,1), 1); hold on;
%     set(b,'FaceColor', [0 0 255]/255, 'EdgeColor', [0 0 0]/255);
%     set(gca, 'FontName', 'Arial')
%     ylabel(xlb{i});
% end
% set(gca, 'FontName', 'Arial')
% xlabel('sample time k');
% 
% printeps(2,'pmpc_lorenz_learnVar'); 

