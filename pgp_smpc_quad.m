function [gppred, syspred, refpred, varargout] = ...
    pgp_smpc_quad(pathtype, mpciterations, systype, varargin)
%%

if nargin > 3 && ~isempty(varargin{1})
    algorithm = varargin{1};
else
    algorithm = 'fminsearch';
end

smpc4quad;

switch algorithm
    case 'fminlbfgs';
        options = setoptimoptions('Algorithm', 'fminlbfgs', 'GradObj','on');
%         options = setoptimoptions('Algorithm', 'fminlbfgs');
    case 'fmincon';
        options = optimoptions('fmincon','Algorithm','active-set',...
            'GradObj','on');   % active-set 
        ubc = umax;
        lbc = umin;
%         lbc = zeros(1, numel(u0))-Inf;
%         ubc = zeros(1, numel(u0))+Inf;
    case 'minimize';
        options = loadoptobj();
    otherwise;
        fprintf('Assigned algorithm is not avaliable. Default fminsearch is used then!\n');
        options = setoptimoptions('Algorithm', 'fminsearch');  
end

tic
fprintf('SMPC Starts...\n');
% main optimization loop
for i=1:N
    ref.x = pathtargets(i:i+H-1,opt.out)';
    ref.u = pathtargets(i:i+H-1,opt.in)';

    costobj = costobjgenerate(dynmodel, ref, H);
    
    switch algorithm
        case 'fminlbfgs';
            [optu, fval(i), exitflag, output] = optiminimize(...
                @(u) MSvalue(u, costobj, m, s), u0, [], [], [], [], [], [], [],...
                options);
        case 'fmincon';
            [optu, fval(i), exitflag, output] = fmincon(...
                @(u) MSvalue(u, costobj, m, s), ...
                u0, [], [], [], [], lbc, ubc, [], options);
        case 'minimize';
            [optu, dfval] = minimize(u0, @MSvalue, options, costobj, m, s);
            fval(i) = dfval(end);
        otherwise;
            [optu, fval(i), exitflag, output] = optiminimize(...
                @(u) MSvalue(u, costobj, m, s), u0, [], [], [], [], [], [], [],...
                options);
    end
    
    u0 = optu;
    [M{i}, S{i}] = pgppred(dynmodel, m, s, u0, H);
    m= M{i}(:,1);
    s = S{i}(:,:,1);
    u = [u, optu(:,1)];
    
    fval(i) = sum((m([1 3 5])-ref.x(1,[1 3 5])').^2);    
    fprintf('SMPC:Iter %6.0d , Fval %12.8d\n', i, fval(i));
end

toc  % count optimization process time

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
        u = normalmatrix('reverse',u',range1(opt.in,:),range2)';
    case 'Gauss01';
        MU = quad_ds.r1(opt.out);
        SIGMA = quad_ds.r2(opt.out);
        res =  repmat(SIGMA, size(res,1),1).*res...
            + repmat(MU, size(res,1),1);
        reff = repmat(SIGMA, size(reff,1),1).*reff...
            + repmat(MU, size(reff,1),1);
        sysres = repmat(SIGMA, size(sysres,1),1).*sysres...
            + repmat(MU, size(sysres,1),1);
        
        MU = quad_ds.r1(opt.in);
        SIGMA = quad_ds.r2(opt.in);
        u = repmat(SIGMA, size(u,1),1).*u...
            + repmat(MU, size(u,1),1);
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

function optparam = loadoptobj()
% parameters 4 minimize optimization process
    optparam.length =  5000;
    optparam.method = 'BFGS';
    optparam.verbosity = 0;
end

function costobj = costobjgenerate(dynmodel, ref, H)
    costobj.dynmodel = dynmodel;
    costobj.H = H;
    costobj.Q =diag([1 1 1 1 1 1]);
    costobj.R =diag([1 1 1]);
    costobj.ref =ref;
end

% f(F, X, varargin{:});
% without gradiences
function J = Mvalue(u, obj, m, s)
    J = 0.0;
    for i =1:obj.H   
        [m, s]= pgppred(obj.dynmodel, m, s, u(:,i), 1);
        J(i) = (m-obj.ref.x(:,i))'*obj.Q*(m-obj.ref.x(:,i));
%         J(i) = (m-obj.ref.x(:,i))'*obj.Q*(m-obj.ref.x(:,i))+trace(obj.Q*s);
    end
    J = sum(J);
end

% f(F, X, varargin{:});
% without consideration of control actions
function [J, dJdu]= MSvalue(u, obj, m, s)
    J = 0.0;
    D = obj.dynmodel.iNum;
    E = obj.dynmodel.sNum;
    for i =1:obj.H  
        [m, s, predobj]= pgppred(obj.dynmodel, m, s, u(:,i), 1);
        
        J(i) = (m-obj.ref.x(:,i))'*obj.Q*(m-obj.ref.x(:,i)) + trace(obj.Q*s);
        dJdm = 2*(m-obj.ref.x(:,i))'*obj.Q;
        dJds = obj.Q;
        
        dmdu = [zeros(E, D); eye(D)];
        dsdu =  zeros((E+D)*(E+D), 3);
        
        dMdu = predobj.dMdm*dmdu + predobj.dMds*dsdu;
        dSdu = predobj.dSdm*dmdu + predobj.dSds*dsdu;

        dJdu(i, :) = dJdm(:)'*dMdu+ dJds(:)'*dSdu;
        
    end
    J = sum(J);
    dJdu = dJdu(:);
end

% f(F, X, varargin{:});
%   gradience + consideration of control actions
function [J, dJdu]= MSUvalue(u, obj, m, s)
    J = 0.0;
    D = obj.dynmodel.iNum;
    E = obj.dynmodel.sNum;
    for i =1:obj.H
        [m, s, predobj]= pgppred(obj.dynmodel, m, s, u(:,i), 1);
        
        J(i) = (m-obj.ref.x(:,i))'*obj.Q*(m-obj.ref.x(:,i))+u(:,i)'*obj.R*u(:,i)...
            +trace(obj.Q*s);
        
        dJdm = 2*(m-obj.ref.x(:,i))'*obj.Q;
        dJds = obj.Q;      
        dJdu = 2*u(:,i)'*obj.R;
        
        dmdu = [zeros(E, D); eye(D)];
        dsdu =  zeros((E+D)*(E+D), 3);
        
        dMdu = predobj.dMdm*dmdu + predobj.dMds*dsdu;
        dSdu = predobj.dSdm*dmdu + predobj.dSds*dsdu;
       
        dJdu(i, :) = dJdm(:)'*dMdu+ dJds(:)'*dSdu + dJdu;
        
    end
    J = sum(J);
    dJdu = dJdu(:);
end

%  m,u are vectors
function [tileM, tileS] = jointDistribution(m, s, u)
    dimS = size(s, 1);
    dimU = size(u,1);
    tileM = [m; u];
    tileS = zeros(dimS+dimU,dimS+dimU);
    tileS(1:dimS, 1:dimS) = s;
end

%  m and s is a scalar,     u may be a matric indicating multiple inputs
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


