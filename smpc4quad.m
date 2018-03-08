defaultLearnPathType = pathtype;

% defining dataset for model learning 
filename4learn1 = strcat('ds_nmpc_', defaultLearnPathType, '_trans');
filename4learn2 = strcat('ds_nmpc_', defaultLearnPathType, '_rotate');

% defining another pathtype to follow
path2follow1 = strcat('ds_nmpc_', pathtype, '_trans');
path2follow2 = strcat('ds_nmpc_', pathtype, '_rotate');

switch systype
    case 'TRAN';
        filename4learn = filename4learn1;
        path2follow = path2follow1;
        loadType = 'original';  % original self0109
        H = 5;
        Ts = 1;
        switch pathtype
            case 'Lorenz';
                umax = [0 2 2];
                umin = [-45 -2 -2]; 
            case 'Elliptical';
                umax = [100 0.2 0.2];
                umin = [0 -0.2 -0.2];  
            otherwise;
                umax = [Inf Inf Inf];
                umin = [-Inf -Inf -Inf];   
        end
    case 'ROTATE';
        filename4learn = filename4learn2;
        path2follow = path2follow2;
        loadType = 'self0109';
        H = 1;
        Ts = 1;
        umax = [0.9 0.9 0.9];
        umin = [0.1 0.1 0.1];
    otherwise;
end

tempRandSequences = randperm(mpciterations);
opt.dsSampleRange = sort(tempRandSequences(1:mpciterations));
opt.dsTestRange = 2:mpciterations;
opt.state = 1:9;
opt.in = 7:9;
opt.out = 10:15;

quad_ds = loadSampleData(filename4learn,loadType,opt);
% numtem = numel(quad_ds);

dynmodel.inputs = quad_ds.proData(opt.dsSampleRange,opt.state);
dynmodel.targets = quad_ds.proData(opt.dsSampleRange,opt.out);
dynmodel.train = @train;
dynmodel. iNum =3;
dynmodel. sNum =6;
dynmodel. oNum =6; 

tic

plant = [];
[dynmodel nlml] = dynmodel.train(dynmodel, plant, 1000);

% 32 -> ascii code of 'blank space'
process_mes = strcat('Successfully Training', 32,  systype, ' GP Model\n');
fprintf(process_mes);
toc   % counting dynamic model learning time

pathobj = loadSampleData(path2follow,loadType,opt);
pathtargets = pathobj.proData;
H=1;

m = pathtargets(1, [1:dynmodel.sNum])';
s= diag([0.1*ones(1, dynmodel.sNum)].^2);
N =mpciterations-10;
u0 = zeros(numel(opt.in), H);
u = u0(:,1);
