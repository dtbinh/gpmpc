
% gp modelling results for quadrotor control
%  translational and rotational subsystems


clear; clc; close all;
set(0, 'DefaultAxesFontName', 'Arial');  
set(0, 'DefaultTextFontName', 'Arial');  

mpciterations = 189;   %  180    199   170  260
loadType = 'self0109';  % self0109    original
trainsize = 189;

tempRandSequences = randperm(mpciterations);
opt.dsSampleRange = sort(tempRandSequences(1:trainsize));
opt.dsTestRange = 1:mpciterations;

opt.state = 1:9;
opt.in = 7:9;
opt.out = 10:15;

pathtype = 'Lorenz';    %  Lorenz  Elliptical
subsys = 'rotate';    %  trans  rotate

switch subsys
    case 'trans';
        filename = strcat('ds_nmpc_', pathtype, '_trans');
        titletxt ={'Translational Subsystem -- Traning Errors',...
            'Traning Size = 189'};
        ylb = {'X[m]','Y[m]','Z[m]'};
    case 'rotate';
        filename = strcat('ds_nmpc_', pathtype, '_rotate');
        titletxt ={'Rotational Subsystem -- Traning Errors',...
            'Traning Size = 189'};
        ylb = {'\phi[rad]','\theta[rad]','\psi[rad]'};
end

trainds = loadSampleData(filename,loadType,opt);
dynmodel.inputs = trainds.proData(opt.dsSampleRange,opt.state);
dynmodel.targets = trainds.proData(opt.dsSampleRange,opt.out);
dynmodel.train = @train;
dynmodel. iNum =3;
dynmodel. sNum =6;
dynmodel. oNum =6; 
tic
plant = [];
[dynmodel nlml] = dynmodel.train(dynmodel, plant, 1000);
toc   % counting dynamic model learning time

traindata.inputs = trainds.proData(opt.dsSampleRange,opt.state);
traindata.targets = trainds.proData(opt.dsSampleRange,opt.out);

testdata.inputs = trainds.proData(opt.dsTestRange,opt.state);
testdata.targets = trainds.proData(opt.dsTestRange,opt.out);

S0size = dynmodel. sNum + dynmodel. iNum;
Ssize = dynmodel. sNum;
s= diag(ones(1,Ssize)*0.1.^2);
s0=  zeros(S0size,S0size,trainsize);
for i=1:trainsize
    s0(1:Ssize,1:Ssize,i)=s;
    [trainpred(i,:) s]= gp0d(dynmodel, traindata.inputs(i,:)', s0(:,:,i));
end

switch loadType
    case 'self0109';
        range1 = trainds.r1;
        range2 = trainds.r2;
        trainpred = normalmatrix('reverse',trainpred,range1(opt.out,:),range2);
        reff = normalmatrix('reverse',traindata.targets,range1(opt.out,:),range2);
    case 'Gauss01';
    otherwise;
        reff = traindata.targets;
end

trainMSE = sum(sum((trainpred-reff).^2))/trainsize

figure(1);
for i=1:3
    subplot(3,1,i);
    plot(trainpred(:,2*i-1)-reff(:,2*i-1), 'k', 'LineWidth', 1.5);
    ylabel(ylb{i});
    xlim([0 trainsize]);
     if i==1
        title(titletxt);
    else
    end
end
xlabel('sample time k');

s= diag([ones(1,Ssize)*0.1].^2);
s1=  zeros(S0size,S0size,mpciterations);
for i=1:mpciterations
    s1(1:Ssize,1:Ssize,i)=s;
    [testpred(i,:) s]= gp0d(dynmodel, testdata.inputs(i,:)', s1(:,:,i));
end

switch loadType
    case 'self0109';
        range1 = trainds.r1;
        range2 = trainds.r2;
        testpred = normalmatrix('reverse',testpred,range1(opt.out,:),range2);
        reff = normalmatrix('reverse',testdata.targets,range1(opt.out,:),range2);
    case 'Gauss01';
    otherwise;
        reff = testdata.targets;
end

testMSE = sum(sum((testpred-reff).^2))/mpciterations

figure(2);
labstr = {' X',' Y',' Z'};
for i=1:3
    subplot(3,1,i);
    j=2*(i-1)+1;
    b1=bar(1:mpciterations, reshape(s1(j,j,:), mpciterations,1), 1); hold on;
    set(b1,'FaceColor', [0 0 255]/255, 'EdgeColor', [0 0 0]/255);
    ylabel(strcat('variances on',labstr{i}));
end
xlabel('sampling time k');
tempvar = sum(reshape(s1(1,1,:), mpciterations,1))+...
    sum(reshape(s1(3,3,:), mpciterations,1))+...
    sum(reshape(s1(5,5,:), mpciterations,1));

avrgVAR = tempvar/189/3

figure(3);
for i=1:3
    subplot(3,1,i);
    plot(testpred(:, 2*i-1), 'b-'); hold on;
    plot(reff(:,2*i-1), 'k--');
end