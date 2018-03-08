clc; clear; close all;
%% Demo of IGP based Quadrotor translational modelling
%
%%
setReleventPath;
loadType = 'original';  % original self0109
% ds_elliptical_trans    ds_singlering_trans ds_helix_trans
filename = 'ds_nmpc_Elliptical';    
tempRandSequences = randperm(200);
opt.dsSampleRange = sort(tempRandSequences(1:100));
opt.dsTestRange = 2:200;
opt.in = [1:9];
opt.on = 10:15;

quad_ds = loadSampleData(filename,loadType,opt);
% numtem = numel(quad_ds);

dynmodel.fcn = @gp2d;  
dynmodel.inputs = quad_ds.proData(opt.dsSampleRange,opt.in);
dynmodel.targets = quad_ds.proData(opt.dsSampleRange,opt.on);
dynmodel.train = @train;

% odei = [1 2 3 4 5 6 7 8 9];                     % varibles for the ode solver
augi = [];                                      % variables to be augmented
dyno = [10 11 12 13 14 15];                   % variables to be predicted (and known to loss)
angi = [];                                    % angle variables
dyni = [1 2 3 4 5 6];          % variables that serve as inputs to the dynamics GP
poli = [1 2 3 4 5 6];          % variables that serve as inputs to the policy
difi = [];            % variables that are learned via differences

plant.noise = diag(ones(1,6)*0.01.^2);            % measurement noise
plant.dt = 0.1;
plant.ctrl = @zoh;                                % controler is zero order hold
% plant.odei = odei;
plant.augi = augi;
plant.angi = angi;
plant.poli = poli;
plant.dyno = dyno;
plant.dyni = dyni;
plant.difi = difi;
plant.prop = @propagated;


mu0 = [0 0 0 0 0 0 ]';                  % initial state mean
s0 = diag([0.1 0.1 0.1 0.1 0.1 0.1].^2);   

% policy.fcn = @(policy,m,s)conCat(@congp,@gSat,policy,m,s);
% policy.maxU = [10 10 10]; 
% [mm ss cc] = gTrig(mu0, s0, plant.angi);          
% mm = [mu0; mm]; cc = s0*cc; ss = [s0 cc; cc' ss]; 
% policy.p.inputs = gaussian(mm(poli), ss(poli,poli), 1)'; 
% policy.p.targets = 0.1*randn(1, length(policy.maxU));    
% policy.p.hyp = ...                                        % initialize policy
%   repmat(log([0.1 0.1 0.1 0.1 0.1 0.1 1 0.01]'), 1,3);    


policy.fcn = @(policy,m,s)conCat(@conlin,@gSat,policy,m,s);
policy.maxU = [10 10 10]; 
policy.p.w = repmat([1 0 1 0 1 0], 3,1);
policy.p.b =[0 0 0]';

opt.length = 100;                        % max. number of line searches
opt.MFEPLS = 30;                         % max. number of function evaluations
                                         % per line search
opt.verbosity = 1;                       % verbosity: specifies how much 
                                         % information is displayed during
                                         % policy learning. Options: 0-3
                                         
cost.fcn = @lossQuad;                       % cost function
diagwei = [1 0 1 0 1 0];
cost.W =  diag(diagwei);
cost.gamma = 1;
                                        

[dynmodel nlml] = dynmodel.train(dynmodel, plant, [300 500]);

H =10;
pathtargets = quad_ds.proData(:,opt.on);
% for i=1:20
%     target = pathtargets(i:i+H-1,:);
%     [policy.p fX3] = minimize(policy.p, 'value', opt, mu0, s0, ...
%         dynmodel, policy, plant, cost, H, target);
%     
%     [M S] = pred(policy, plant, dynmodel, mu0, s0, 1);
%     mu0 = M(:, end);
%     s0 = S(:,:,end);
%     
%     outobj{i}.policy = policy;
%     outobj{i}.m = mu0;
%     outobj{i}.s0 = s0;
% end

[M S] = pred(policy, plant, dynmodel, mu0, s0, 200);

% opt.fh = 1;
% [policy.p fX3] = minimize(policy.p, 'value', opt, mu0, s0, ...
%   dynmodel, policy, plant, cost, H);

y_pre = M'; y_pre = y_pre(2:end, :);
y_test = quad_ds.proData(:,opt.on);

% figure(1);
% plot(M(1,:), 'k-'); hold on;
% plot(M(2,:), 'b*-'); hold on;
% plot(M(3,:), 'go-'); hold on;
% plot(M(4,:), 'r^-'); hold on;
% plot(M(5,:), 'b:'); hold on;
% plot(M(6,:), 'y-'); hold on;

% switch loadType
%     case 'self0109';
%         range1 = quad_ds.r1;
%         range2 = quad_ds.r2;
%         y_pre = normalmatrix('reverse',y_pre,range1(opt.on,:),range2);
%         y_test = normalmatrix('reverse',y_test,range1(opt.on,:),range2);
%     case 'Gauss01';
%         y_test = y_test;
% end
% 
% yplot = [y_pre(:,1), y_pre(:,3), y_pre(:,5)];
% yrplot = [y_test(:,1), y_test(:,3), y_test(:,5)];
% 
% plotopt.xlabel = 'Times';
% plotopt.ylabel = {'x','y','z'};
% plotopt.yerrlabel = {'','',''};
% plotopt.realcurvetype = 'b--';
% plotopt.predcurvetype = 'k-';
% plotopt.errcurvetype ='k--';
% plotopt.reallegend ='Translational Function';
% plotopt.predlegend = 'CGP model';
% plotopt.errlegend = 'err';
% 
% % gpplot(yrplot, yplot, var);
% figure(1);
% drawRealPrediction(yrplot,yplot,0,plotopt);


