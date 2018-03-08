function [gppred, refpred, varargout] = ...
    pgp_smpc(dynmodel, mpcobj, systype)
%%

N =mpcobj.N;
H = mpcobj.H;
pathtargets = dynmodel.ds.proData;
opt = dynmodel.ds.opt;

plant = [];
[dynmodel nlml] = dynmodel.train(dynmodel, plant, [300 500]);

% 32 -> ascii code of 'blank space'
process_mes = strcat('Successfully Training', 32,  systype, ' GP Model\n');
fprintf(process_mes);

m = zeros(dynmodel.iNum+dynmodel.sNum,1);
s = diag([0.01+zeros(1, dynmodel.sNum), zeros(1,dynmodel.iNum)]);
u = [];
u0 = zeros(numel(opt.in), H);
lbc = zeros(1, numel(u0))-Inf;
ubc = zeros(1, numel(u0))+Inf;

fprintf('SMPC Starts...\n');
for i=1:N
    ref.x = pathtargets(i:i+H-1,opt.out)';
    ref.u = pathtargets(i:i+H-1,opt.in)';
%     [M{i} S{i}] = pgppred(dynmodel, m(:,i), s(:,:,i), H);
    costobj = costobjgenerate(dynmodel, ref, H);
    [optu, fval(i), exitflag, output] = fmincon(...
        @(u) probabilisticost(costobj, m(:,i), s(:,:,i), u), ...
        u0, [], [], [], [], lbc, ubc);
    u = [u, optu(:,1)];
    u0 = optu;
    m(dynmodel.sNum+1:end, i) = u(:,end);
    [M{i}, S{i}] = pgppred(dynmodel, m(:,i), s(:,:,i), 1);
    m(1:dynmodel.sNum, i+1) = M{i}(:, end);
    s(1:dynmodel.sNum, 1:dynmodel.sNum, i+1) = S{i}(:, :, end);
    fprintf('SMPC:Iter %6.0d , Fval %12.8d\n', i, fval(i));
end

for j=1:N
    res(j,:) = pgppred(dynmodel, m(:,j), s(:,:,j), 1);    
end

reff = pathtargets(:,opt.out);

switch dynmodel.loadType
    case 'self0109';
        range1 = dynmodel.ds.r1;
        range2 = dynmodel.ds.r2;
        res = normalmatrix('reverse',res,range1(opt.out,:),range2);
        reff = normalmatrix('reverse',reff,range1(opt.out,:),range2);
    case 'Gauss01';
    otherwise;
end

gppred = res;
refpred = reff;

if nargout >2
    varargout{1} = u;
    varargout{2} = M;
    varargout{3} = S;
    varargout{4} = fval;
end
end

function costobj = costobjgenerate(dynmodel, ref, H)
    costobj.dynmodel = dynmodel;
    costobj.H = H;
    costobj.Q =diag(1+zeros(1, dynmodel.oNum));
    costobj.R =diag(1+zeros(1, dynmodel.iNum));
    costobj.ref =ref;
end

function cost = probabilisticost(obj, m, s, u)
    Q = obj.Q;   R = obj.R;  ref = obj.ref; 
    H = obj.H;
    dyn = obj.dynmodel;
    
     cost = 0.0;
    for i =1:H
        m(dyn.sNum+1:end) = u(:,i);
        [mu(:,i) var(:,:,i)]= pgppred(dyn, m, s, 1);
        m(1:dyn.sNum) = mu(:,i);
        s(1:dyn.sNum, 1:dyn.sNum)=var(:,:,i);
        cost = cost + (mu(:,i)-ref.x(:,i))'*Q*(mu(:,i)-ref.x(:,i)) +...
            u(:,i)'*R*u(:,i);
    end
    cost = cost + trace(Q*var(:,:,end));
end

function [M, S] = pgppred(dynmodel, m, s, H)
    for i=1:H
        sNum = dynmodel.sNum;
        [m(1:sNum) s(1:sNum, 1:sNum)] = gp2d(dynmodel, m, s);
        M(:,i) = m(1:sNum) ;
        S(:,:,i) = s(1:sNum, 1:sNum);
    end
end

