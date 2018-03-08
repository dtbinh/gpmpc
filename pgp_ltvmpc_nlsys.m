function  [gppred, varargout] = pgp_ltvmpc_nlsys(varargin)
%
%   ds=[x1,x2,x3,x4, u1, u2, x1,x2,x3,x4]
%   y1 =x1; y2 = x3;
%
%%  step control result is not good.  (no idea why)
%  For ltv system control,  fminsearch performs better than minimize  
%%
if nargin >0 && ~isempty(varargin{1})
    syscase = varargin{1};
else
    syscase = 'nltv';
end

if nargin >1 && ~isempty(varargin{2})
    optimizer = varargin{2};
else
    optimizer = 'ltvmpc';
end

if nargin >2 && ~isempty(varargin{3})
    loadType = varargin{3};
else
    loadType = 'original';
end

if nargin>3 && ~isempty(varargin{4})
    pathcase = varargin{4};
else
    pathcase = 'lorenz';
end

switch syscase
    case 'ltv';
        dynmodel. iNum =2;
        dynmodel. sNum =2;
        dynmodel. oNum =2; 
        opt.state = 1:4;
        opt.in = 3:4;
        opt.out = 5:6;
        switch pathcase
            case 'duffing';
                filename = 'lmpc4ltvsys';
                mpciterations = 250;
                umax = [Inf Inf];
                umin = [-Inf -Inf];
            case 'lorenz';
                filename = 'lmpc4ltv_lorenz';
                mpciterations = 180;
                umax = [Inf Inf];
                umin = [-Inf -Inf];
        end
    case 'nltv';
        dynmodel. iNum =2;
        dynmodel. sNum =4;
        dynmodel. oNum =4; 
        opt.state = 1:6;
        opt.in = 5:6;
        opt.out = 7:10;
        opt.outindex = [7 9];
%         pathcase = 'duffing';  %  duffing  lorenz  sincos
        switch pathcase
            case 'duffing';
                filename = 'ds_nmpc_nltv_duffing';
                mpciterations = 199;
                umax = [0.25 0.2];
                umin = [-0.5 -0.5];
            case 'lorenz';
                filename = 'ds_nmpc_nltv_lorenz';
                mpciterations = 189; 
                umax = [4 7];
                umin = [-4 -7];
            case 'sincos';
                filename = 'ds_nmpc_nltv_sincos';
                mpciterations = 170;
                umax = [1.5 1.5];
                umin = [-1.5 -1.5];
            case 'step';
                filename = 'ds_nmpc_nltv_step';
                mpciterations = 199;
                umax = [1.75 2.5];
                umin = [0.25 0.5];
        end
end

% umax = [Inf Inf];
% umin = [-Inf -Inf];

tempRandSequences = randperm(mpciterations);
opt.dsSampleRange = sort(tempRandSequences(1:round(mpciterations*1)));
opt.dsTestRange = 2:mpciterations;

trainds = loadSampleData(filename,loadType,opt);
dynmodel.inputs = trainds.proData(opt.dsSampleRange,opt.state);
dynmodel.targets = trainds.proData(opt.dsSampleRange,opt.out);
dynmodel.train = @train;
tic

plant = [];
[dynmodel nlml] = dynmodel.train(dynmodel, plant, 1000);

% 32 -> ascii code of 'blank space'
process_mes = strcat('Successfully Training', 32, ' GP Model\n');
fprintf(process_mes);
toc   % counting dynamic model learning time

H=10;
pathtargets = trainds.proData;
m = trainds.proData(1, 1:dynmodel.sNum)';
s= diag([0.1*ones(1, dynmodel.sNum)].^2);
N =mpciterations-10;
u = zeros(dynmodel.iNum, 1);
deltaU = zeros(dynmodel.iNum*H,1);
switch optimizer
    case 'fminlbfgs';
        options = setoptimoptions('Algorithm', 'fminlbfgs', 'GradObj','on');
        options = setoptimoptions('Algorithm', 'fminlbfgs');
        H =1;
%         options = setoptimoptions('Algorithm', 'fminlbfgs');
    case 'fminsearch';
        options = setoptimoptions('Algorithm', 'fminsearch');
        H =1;
    case 'fmincon';
        options = optimoptions('fmincon','Algorithm','sqp');   % active-set 
        H =1;
        lbc = umin;
        ubc = umax;
    case 'minimize';
        options = loadoptobj();
        H =1;
%     otherwise;
%         fprintf('Assigned algorithm is not avaliable. Default fminsearch is used then!\n');
%         options = setoptimoptions('Algorithm', 'fminsearch'); 
    otherwise;
end

tic
fprintf('GPMPC Starts...\n');
% main optimization loop
for i=1:N
    ref.x = pathtargets(i:i+H-1,opt.out)';
    ref.u = pathtargets(i:i+H-1,opt.in)';

    costobj = costobjgenerate(dynmodel, ref, H);
    
    switch optimizer
        case {'fminlbfgs', 'fminsearch'};    % H=1  only
            [optu, fval(i), exitflag, output] = optiminimize(...
                @(u) MSUvalue(u, costobj, m, s), u(:,i), [], [], [], [], [], [], [],...
                options);
        case 'minimize';   % H=1  only
            [optu, dfval] = minimize(u(:,i), @MSUvalue, options, costobj, m, s);
            fval(i) = dfval(end);   
        case 'fmincon';   % H=1  only
            [optu, fval(i), exitflag, output] = fmincon(...
                @(u) MSUvalue(u, costobj, m, s), ...
                u(:,i), [], [], [], [], lbc, ubc, [], options);
        case 'quadprog';
            ref.x = ref.x';
            ref.u = ref.u';
            [m0, s0, predobj]= pgppred(dynmodel, m, s, u(:,i), 1);
%----------------------------------------------------------------------
%   xk = muk
%             [ltvm.A, ltvm.B]= linearModel(dynmodel, predobj.dMdm, predobj.dMds);  
%             Q = diag([1 1 1 1]);
%----------------------------------------------------------------------
%----------------------------------------------------------------------
%   xk = [muk, vec(sigk)]
            [ltvm.A, ltvm.B]= linearProbmodel(dynmodel, predobj.dMdm, ...
                predobj.dMds, predobj.dSdm, predobj.dSds,m0,s0);  
            [dA, dB] = c2d(ltvm.A, ltvm.B, 0.01);
%             dA = ltvm.A;
%             dB = ltvm.B;
            Q = diag(ones(dynmodel.sNum+dynmodel.sNum^2,1));
%----------------------------------------------------------------------
            R = diag([27000 21000]);
            [tildeA, tildeB, tildeQ, tildeR] =...
                extendMatrix(H, dA, dB, Q, R);
            Bk = tildeB'*tildeQ*tildeB+ tildeR;
%             [HBk,exitflag] = chol(Bk);
%             if(exitflag)
%                 [Hes, jitter] = jit_chol(Bk);
%             else
%                 Hes = HBk;
%             end

            Hes = Bk;
            setp = [ref.x, repmat(0*s(:)', H, 1)];
            x0 =[m; s(:)]-setp(1,:)';
            f = 2*(x0)'*tildeA'*tildeQ*tildeB;
            lbc = [umin(1).*ones(H,1)-ref.u(:,1),...
                umin(2).*ones(H,1)-ref.u(:,2)];
            ubc = [umax(1).*ones(H,1)-ref.u(:,1),...
                umax(2).*ones(H,1)-ref.u(:,2)];
            lbc = reshape(lbc', H*dynmodel.iNum,1);
            ubc = reshape(ubc', H*dynmodel.iNum,1);
            
            optuu= quadprog(Hes,f,[],[],[],[],lbc,ubc,...
                deltaU-reshape(ref.u', H*dynmodel.iNum,1));
            optu= reshape(optuu', dynmodel.iNum,H)';
            optu = optu(1,:);
        case 'sqp';
            ref.x = ref.x';
            ref.u = ref.u';
            [m0, s0, predobj]= pgppred(dynmodel, m, s, u(:,i), 1);
            [ltvm.A, ltvm.B]= linearModel(dynmodel, predobj.dMdm, predobj.dMds);  
            [dA, dB] = c2d(ltvm.A, ltvm.B,  0.01);
            Q = diag([1 1 1 1]);
            R = diag([27000 21000]);
            [tildeA, tildeB, tildeQ, tildeR] =...
                extendMatrix(H, dA,dB, Q, R);
            Bk = tildeB'*tildeQ*tildeB+ tildeR;
            
            f = 2*(m-ref.x(1,:)')'*tildeA'*tildeQ*tildeB;
%             Aeq = ltvm.B;
%             Beq = m0-temprefx-ltvm.A*(m-ref.x');
%             optuu= quadprog(Hes,f,[],[],[],[],lbc,ubc,u(:,i)-ref.u');
%             optuu= quadprog(Hes,f,[],[],Aeq,Beq,lbc,ubc,u(:,i)-ref.u');
%             uu0 = u(:,i)-reshape(ref.u', H*dynmodel.iNum, 1);
%             lbc = [umin(1)-ref.u(1,1) umin(2)-ref.u(1,2)];
%             ubc = [umax(1)-ref.u(1,1)  umax(2)-ref.u(1,2)];
%             lbc = repmat(lbc, H,1);
%             ubc = repmat(ubc, H,1);
            lbc = [umin(1)*ones(H,1)-ref.u(:,1),...
                umin(2)*ones(H,1)-ref.u(:,2)];
            ubc = [umax(1)*ones(H,1)-ref.u(:,1),...
                umax(2)*ones(H,1)-ref.u(:,2)];
            
            lbc = reshape(lbc', H*dynmodel.iNum,1);
            ubc = reshape(ubc', H*dynmodel.iNum,1);
            options = optimoptions('fmincon','Algorithm','sqp'); 
            [optuu, fval(i)] = fmincon(...
                @(u) sqpvalue(Bk,f,u,s0), ...
                deltaU-reshape(ref.u', H*dynmodel.iNum,1), [], [], [], [], lbc,ubc, [],options);
            optu= reshape(optuu', dynmodel.iNum,H)';
            optu = optu(1,:);
        otherwise;
            % using extend gp model based linear mpc
            ref.x = ref.x';
            ref.u = ref.u';
            init0 =[m; s(:)];
            [m0, s0, predobj]= pgppred(dynmodel, m, s, u(:,i), 1);
            [ltvm.A, ltvm.B]= linearProbmodel(dynmodel, predobj.dMdm, ...
                predobj.dMds, predobj.dSdm, predobj.dSds,m0,s0);  
            ltvm.C = diag(ones(length(m)+length(s(:)), 1));  
            ltvm.D=0;
            
            con.u = [umin(1)-ref.u(1,1), umax(1)-ref.u(1,1), Inf;...
                umin(2)-ref.u(1,2), umax(2)-ref.u(1,2), Inf];
%             con.u = [];
            con.y = [];

            Ts=0.1;
            plant = jSS(ltvm.A, ltvm.B, ltvm.C, ltvm.D);
            plant = c2d(plant, Ts);

            model = plant;
            plant.x0 =  init0;
            model.x0 =  init0;

            uwt = [21000 27000]';
            ywt = ones(length(m)+length(s(:)), 1);
           
            setp = [ref.x, repmat(0*s(:)', H, 1)];
            MPC = jMPC(model, H, H, uwt,ywt,con);
 
            simopts = jSIM(MPC,plant, H, setp);
            simresult = sim(MPC,simopts,'Mex');
            optu =  simresult.plotvec.u(2,:);
            yp =  simresult.plotvec.yp(2,:);
            jmpcU(i,:) = optu;
    end

    switch optimizer
        case {'fminlbfgs', 'fminsearch', 'minimize', 'fmincon'};
            u = [u, optu];
            [m, s]= pgppred(dynmodel, m, s, u(:,i+1), 1);
            Mmat(:,i) = m;
            Smat(:,:,i) =s;
%             jmpcU=0;
        otherwise;
            u = [u, optu'+ref.u(1,:)'];
            [m, s]= pgppred(dynmodel, m, s, u(1:dynmodel.iNum,i+1), 1);
            Mmat(:,i) = m;
            Smat(:,:,i) =s;
            switch syscase
                case 'ltv';
                    fval(i) = sum((m-ref.x(1,:)').^2);   % for ltv
                case 'nltv';
                    fval(i) = sum((m([1 3])-ref.x(1,[1 3])').^2);  % for nltv
            end
    end
    fprintf('GPMPC:Iter %6.0d , Fval %12.8d\n', i, fval(i));
end

toc  % count optimization process time

res = Mmat';
reff = pathtargets(:,opt.out);

switch loadType
    case 'self0109';
        range1 = trainds.r1;
        range2 = trainds.r2;
        res = normalmatrix('reverse',res,range1(opt.out,:),range2);
        reff = normalmatrix('reverse',reff,range1(opt.out,:),range2);
    case 'Gauss01';
    otherwise;
end

gppred = res;
refpred = reff;
% mobj.pathtargets = pathtargets;
% mobj.jmpcU = jmpcU;

if nargout >1
    varargout{1} = refpred;
    varargout{2} = u;
    varargout{3} = fval;
    varargout{4} = Smat;
%     varargout{5} = mobj;
end

end

function [J, dJdu] = sqpvalue(H,f,u,s)
    J = u'*H*u+f*u + 10*trace(s);
    dJdu = 0;
end

function optparam = loadoptobj()
% parameters 4 minimize optimization process
    optparam.length =  1000;
    optparam.method = 'BFGS';
    optparam.verbosity = 0;
end

function costobj = costobjgenerate(dynmodel, ref, H)
    costobj.dynmodel = dynmodel;
    costobj.H = H;
    costobj.Q =diag(ones(dynmodel. sNum,1));
    costobj.R =diag(ones(dynmodel. iNum,1));
    costobj.ref =ref;
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
        dsdu =  zeros((E+D)*(E+D), D);
        
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
        dsdu =  zeros((E+D)*(E+D), D);
        
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

function [A, B] = ...
    linearProbmodel(dynmodel, dMdm, dMds, dSdm, dSds, m0, s0)
    D = dynmodel.iNum;
    E = dynmodel.sNum;
    
    dS0dsqrtS0 = 2*repmat(reshape(s0^(0.5), 1, E^2)', 1,E^2);
    dS0dsqrtS0 = 1;
    dS0S0 = 2*s0^(0.5);
    
    dmdM0 =  [eye(E); zeros(D, E)];
    dsdS0 = [eye(E*E); zeros(D*D+2*E*D, E*E)];
    dMdM0 = dMdm*dmdM0; 
    dMdS0 = dMds*dsdS0;       dMdS0 = dMdS0*dS0dsqrtS0;
    dSdM0 = dSdm*dmdM0;    dSdM0 =  (1./dS0dsqrtS0)'*dSdM0;
    dSdS0 = dSds*dsdS0;
    
    A = [dMdM0, dMdS0; dSdM0, dSdS0];
    
    dmdu = [zeros(E, D); eye(D)];
    dsdu = zeros((E+D)*(E+D), D);
    dMdu = dMdm*dmdu + dMds*dsdu;
    dSdu = dSdm*dmdu + dSds*dsdu;  dSdu =  (1./dS0dsqrtS0)'*dSdu;
    
    B = [dMdu; dSdu];
end

function [A, B] = ...
    linearModel(dynmodel, dMdm, dMds)
    D = dynmodel.iNum;
    E = dynmodel.sNum;
    
    dmdM0 =  [eye(E); zeros(D, E)];
    dMdM0 = dMdm*dmdM0; 
    A = dMdM0;
    
    dmdu = [zeros(E, D); eye(D)];
%     dsdu = zeros((E+D)*(E+D), D);
%     dMdu = dMdm*dmdu + dMds*dsdu;
    dMdu = dMdm*dmdu;
    
    B = dMdu;
end

function [tildeA, tildeB, tildeQ, tildeR] = extendMatrix(H, A, B, Q, R)
    tildeA = cell(H,1);
    for i =1:H
         tildeA{i} = A^H;
    end
    tildeA = cell2mat(tildeA);
    tildeB = cell(H,H);
    for j = 1:H
        for i= 1:H
            if i>=j
                tildeB{i,j}= A^(i-j)*B;
            else
                tildeB{i,j}= 0.*A*B;
            end
            
        end
    end
    tildeB = cell2mat(tildeB);
    tildeB = tril(tildeB);
    
    tildeQ = diag(repmat(diag(Q),H,1));
    tildeR = diag(repmat(diag(R),H,1));
end
