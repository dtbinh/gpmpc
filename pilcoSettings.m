tempRandSequences = randperm(mpciterations);
opt.dsSampleRange = sort(tempRandSequences(1:mpciterations));
opt.dsTestRange = 2:mpciterations;
opt.state = 1:9;
opt.in = 7:9;
opt.out = 10:15;
quad_ds = loadSampleData(filename,loadType,opt);

odei = [1 2 3 4 5 6 7 8 9];            % varibles for the ode solver
augi = [];                   % variables to be augmented
dyno = [10 11 12 13 14 15];            % variables to be predicted (and known to loss)
angi = [];                  % angle variables
dyni = [1 2 3 4 5 6];          % variables that serve as inputs to the dynamics GP
poli = [1 2 3 4 5 6];          % variables that serve as inputs to the policy
difi = [];            % variables that are learned via differences

plant.noise = diag(ones(1,6)*0.1.^2);            % measurement noise
plant.dt = 0.10;
plant.ctrl = @zoh;                                % controler is zero order hold
plant.odei = odei;
plant.augi = augi;
plant.angi = angi;
plant.poli = poli;
plant.dyno = dyno;
plant.dyni = dyni;
plant.difi = difi;
plant.prop = @propagated;

policy.fcn = @conlin;% controller 
policy.maxU = [100 100 100];      % max. amplitude of 
policy.p.w = repmat([1 1 1 1 1 1], 3,1);   %  u= w*X+b;
policy.p.b = [0 0 0]';

% 5. Set up the cost structure
cost.fcn = @lossQuad;                  % cost function
cost.gamma = 1;                            % discount factor
cost.W = diag([1 1 1 1 1 1].^2);

dynmodel.fcn = @gp0d;                % function for GP predictions
dynmodel.train = @train;             % function to train dynamics model
dynmodel.inputs = quad_ds.proData(opt.dsSampleRange,opt.state);
dynmodel.targets = quad_ds.proData(opt.dsSampleRange,opt.out);

optparam.length =  100;
optparam.method = 'LBFGS';
optparam.verbosity = 1;
% optparam.mem = 100;