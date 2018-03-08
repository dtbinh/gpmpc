
clear; clc; close all;
setReleventPath;

load samples.txt;
dsx = samples(:,1:8);
dsy = samples(:,9:10);

treatType= 3;
switch treatType
    case 1;
    case 2;
    case 3;
        [ddsx, Psx]=mapminmax(dsx', 0.1,0.9); 
        [ddsy, Psy]=mapminmax(dsy',  0.1,0.9);
        ddsx= ddsx';
        ddsy= ddsy';
    otherwise;
        ddsx = dsx;
        ddsy = dsy;
end

x_test=ddsx(41+1:end,:);
y_test=ddsy(41+1:end,:);

dynmodel.fcn = @gp0;  
dynmodel.inputs = ddsx(1:41,:);
dynmodel.targets = ddsy(1:41,:);
dynmodel.train = @train;

S0 = diag([0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01].^2);   
plant = [];
[dynmodel nlml] = dynmodel.train(dynmodel, plant, [300 500]);

for i=1:size(x_test, 1)
    M{i} = gp2d(dynmodel, x_test(i,:)', S0);
end

y_pre = cell2mat(M)';

switch treatType
    case 1;
    case 2;
    case 3;
        revpre = mapminmax('reverse', y_pre',Psy);
        revtest = mapminmax('reverse',y_test',Psy);
    otherwise;
        revpre = y_pre';
        revtest = y_test';
end

plotytest = revtest';
plotypre = revpre';
figure(1); 
sset.range= {[300,600],[400,700], [10,40]};
sset.label ={'YS','TS','EL'};
sset.lineinterval = {2,2,0.4};
for i=1:size(dsy,2)
    if size(dsy,2)<=2
        subplot(1,size(dsy,2),i);
    elseif size(dsy,2)<=4
        subplot(2,2,i);
    end
    DataFitCurve(plotytest(:,i),plotypre(:,i),...
        sset.lineinterval{i}, sset.range{i},sset.label{i}); hold on;
end

figure(2);
res1 = showErrorStatistic(plotytest,plotypre,'RMSE');
res2 = showErrorStatistic(plotytest,plotypre,'R2');
for i=1:size(dsy,2)
    subplot(size(dsy,2),1,i);
    plot(plotytest(:,i),'k*-'); hold on;
    plot(plotypre(:,i),'bo-'); 
    titlestr1 = strcat('RMSE=',num2str(res1{i}));
    titlestr2 = strcat(' R^2=',num2str(res2{i}));
    title({titlestr1,titlestr2},'FontSize', 8, 'FontWeight', 'Bold');
end


