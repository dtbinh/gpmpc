function [gppred, syspred, refpred, varargout] = ...
    pgp_lqrsmpc_quad(pathtype, mpciterations, systype, varargin)
%%
%  stochastic mpc control using lqr based on local linearization technique
%%

if nargin > 3 && ~isempty(varargin{1})
    algorithm = varargin{1};
else
    algorithm = 'fminsearch';
end

filename1 = strcat('ds_nmpc_', pathtype, '_trans');
filename2 = strcat('ds_nmpc_', pathtype, '_rotate');

switch systype
    case 'TRAN';
        filename = filename1;
        loadType = 'original';  % original self0109
        H = 1; 
    case 'ROTATE';
        filename = filename2;
        loadType = 'self0109';
        H = 1;
    otherwise;
end

tempRandSequences = randperm(mpciterations);
opt.dsSampleRange = sort(tempRandSequences(1:mpciterations));
opt.dsTestRange = 2:mpciterations;
opt.state = 1:9;
opt.in = 7:9;
opt.out = 10:15;

quad_ds = loadSampleData(filename,loadType,opt);
pathtargets = quad_ds.proData;
% numtem = numel(quad_ds);

dynmodel.inputs = quad_ds.proData(opt.dsSampleRange,opt.state);
dynmodel.targets = quad_ds.proData(opt.dsSampleRange,opt.out);
dynmodel.train = @train;
dynmodel. iNum =3;
dynmodel. sNum =6;
dynmodel. oNum =6; 

tic
[dynmodel nlml] = dynmodel.train(dynmodel, [], 500);
% 32 -> ascii code of 'blank space'
process_mes = strcat('Successfully Training', 32,  systype, ' GP Model\n');
fprintf(process_mes);
toc   % counting dynamic model learning time

m = pathtargets(1, [1:dynmodel.sNum])';
s = diag([0.1 0.1 0.1 0.1 0.1 0.1].^2);
N =mpciterations-10;
u0 = zeros(numel(opt.in), H);
u = u0(:,1);
tic

lbc = zeros(1, numel(u0))-Inf;
ubc = zeros(1, numel(u0))+Inf;
options = setoptimoptions('Algorithm', 'fminsearch');

fprintf('SMPC Starts...\n');
for i=1:N
    ref.x = pathtargets(i:i+H-1,opt.out)';
    ref.u = pathtargets(i:i+H-1,opt.in)';
    costobj = costobjgenerate(dynmodel, ref, H);
   
    [optu, fval(i), exitflag, output] = optiminimize(...
        @(u) MSvalue(u, costobj, m, s), u0, [], [], [], [], [], [], [],...
        options);

    [M{i}, S{i}] = pgppred(dynmodel, m, s, u0, H);
    m= M{i}(:,1);
    s = S{i}(:,:,1);
    u = [u, optu(:,1)];
    u0 = optu;
    
    fprintf('SMPC:Iter %6.0d , Fval %12.8d\n', i, fval(i));
end
toc

GNobj = GNgenerate(200+5, 3, 0, 0.01);

for j=1:N
    switch systype
        case 'TRAN';
            sysres(j,:) = subsys_trans(M{j}(:,1), u(:, j), GNobj.GNnoise(j, :));
            M{j} = M{j}(:, 1);
        case 'ROTATE';
            sysres(j,:) = subsys_rotate(M{j}(:,1), u(:, j), GNobj.GNnoise(j, :));
            M{j} = M{j}(:, 1);
        otherwise;
    end
end

res = cell2mat(M)';
reff = pathtargets(:,opt.out);

switch loadType
    case 'self0109';
        range1 = quad_ds.r1;
        range2 = quad_ds.r2;
        res = normalmatrix('reverse',res,range1(opt.out,:),range2);
        sysres = normalmatrix('reverse',sysres,range1(opt.out,:),range2);
        reff = normalmatrix('reverse',reff,range1(opt.out,:),range2);
    case 'Gauss01';
    otherwise;
end

gppred = res(:,[1 3 5]);
syspred = sysres(:,[1 3 5]);
refpred = reff(:,[1 3 5]);

if nargout >3
    varargout{1} = u;
    varargout{2} = M;
    varargout{3} = S;
    varargout{4} = fval;
end
end

function costobj = costobjgenerate(dynmodel, ref, H)
    costobj.dynmodel = dynmodel;
    costobj.H = H;
    costobj.Q =diag([1 1 1 1 1 1]);
    costobj.R =diag([1 1 1]);
    costobj.ref =ref;
end

function [A, B] = ...
    linearProbmodel(dynmodel, dMdm, dMds, dSdm, dSds)
    D = dynmodel.iNum;
    E = dynmodel.sNum;
    
    dmdM0 =  [eye(E); zeros(D, E)];
    dsdS0 = [eye(E*E); zeros(D*D+2*E*D, E*E)];
    dMdM0 = dMdm*dmdM0;
    dMdS0 = dMds*dsdS0;
    dSdM0 = dSdm*dmdM0;
    dSdS0 = dSds*dsdS0;
    
    A = [dMdM0, dMdS0; dSdM0, dSdS0];
    
    dmdu = [zeros(E, D); eye(D)];
    dsdu = zeros((E+D)*(E+D), 3);
    dMdu = dMdm*dmdu + dMds*dsdu;
    dSdu = dSdm*dmdu + dSds*dsdu;
    
    B = [dMdu; dSdu];
end

function [M, S, varargout] =...
    pgppred(dynmodel, m, s, u, H)

    for i=1:H
        [tileM, tileS] = jointDistribution(m, s, u(:,i));
        [m s, V, dMdm, dSdm, dVdm, dMds, dSds, dVds] = gp0d(dynmodel, tileM, tileS);
        M(:,i) = m;
        S(:,:,i) = s;
        predobj{i}.V = V;
        predobj{i}.dMdm = dMdm;
        predobj{i}.dSdm = dSdm;
        predobj{i}.dVdm = dVdm;
        predobj{i}.dMds = dMds;
        predobj{i}.dSds = dSds;
        predobj{i}.dVds = dVds;
    end
    
    if nargout > 2
        varargout = predobj;
    end
end

function J= MSvalue(u, obj, m, s)
    J = 0.0;
    for i =1:obj.H  
        [m, s, predobj]= pgppred(obj.dynmodel, m, s, u(:,i), 1);
        J(i) = (m-obj.ref.x(:,i))'*obj.Q*(m-obj.ref.x(:,i)) + trace(obj.Q*s);
    end
    [A, B]= linearProbmodel(obj.dynmodel, predobj.dMdm, ...
        predobj.dMds, predobj.dSdm, predobj.dSds);
    Qu = diag(ones(size(A, 2), 1));
    Ru = diag(ones(size(B, 2), 1));
    [Ku, Pu, Eu] = dlqr(A, B, Qu, Ru, 0);
    mm = [m; s(:)];
    J = sum(J) + mm'*Pu*mm;
end

function [tileM, tileS] = jointDistribution(m, s, u)
    dimS = size(s, 1);
    dimU = size(u,1);
    tileM = [m; u];
    tileS = zeros(dimS+dimU,dimS+dimU);
    tileS(1:dimS, 1:dimS) = s;
end


function y = subsys_trans(x, u, noise)
    
    phi = pi/4;  theta = pi/4;
    m = 0.74;    g = 9.81;
    y = zeros(6, 1);

    y(1) = x(2);
    y(2) = u(2)*u(1)/m + noise(1);
    y(3) = x(4);
    y(4) = u(3)*u(1)/m+ noise(2);
    y(5) = x(6);
    y(6) = g-(cos(phi)*cos(theta))*u(1)/m+ noise(3);
end

function y = subsys_rotate(x, u, noise)

    quad_params;
    y = zeros(6, 1);
    y(1) = x(2);
    y(2) = x(4)*x(6)*a1+x(4)*a2*OmegaR+b1*u(1)+ noise(1);
    y(3) = x(4);
    y(4) = x(2)*x(6)*a3-x(2)*a4*OmegaR+b2*u(2)+ noise(2);
    y(5) = x(6);
    y(6) = x(4)*x(2)*a5+b3*u(3)+ noise(3);
end