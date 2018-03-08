function [gppred, varargout] = ...
    pgp_ltvmpc_quad(pathtype, mpciterations, systype)
%%
%  stochastic mpc control using lqr based on local linearization technique
%%

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
        H = 10;
        Ts = 0.01;
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
        H = 10;
        Ts = 1;
        umax = [0.9 0.9 0.9];
        umin = [0.1 0.1 0.1];
    otherwise;
end

tempRandSequences = randperm(mpciterations);
opt.dsSampleRange = sort(tempRandSequences(1:100));
opt.dsTestRange = 2:mpciterations;
opt.state = 1:9;
opt.in = 7:9;
opt.out = 10:15;

quad_ds = loadSampleData(filename4learn,loadType,opt);
% pathtargets = quad_ds.proData;
% numtem = numel(quad_ds);

dynmodel.inputs = quad_ds.proData(opt.dsSampleRange,opt.state);
dynmodel.targets = quad_ds.proData(opt.dsSampleRange,opt.out);
dynmodel.train = @train;
dynmodel. iNum =3;
dynmodel. sNum =6;
dynmodel. oNum =6; 

tic
[dynmodel nlml] = dynmodel.train(dynmodel, [], 1000);
% 32 -> ascii code of 'blank space'
process_mes = strcat('Successfully Training', 32,  systype, ' GP Model\n');
fprintf(process_mes);
toc   % counting dynamic model learning time


pathobj = loadSampleData(path2follow,loadType,opt);
pathtargets = pathobj.proData;

controller = 'mycontroller';    %  mycontroller     jMPC

m = pathtargets(1, 1:dynmodel.sNum)';
s= diag([0.1*ones(1, dynmodel.sNum)].^2);
N =mpciterations-10;
u = zeros(dynmodel.iNum, 1);
deltaU = zeros(dynmodel.iNum*H,1);
tic

fprintf('LTVMPC Starts...\n');
for i=1:N
    ref.x = pathtargets(i:i+H-1,opt.out);
    ref.u = pathtargets(i:i+H-1,opt.in);
      
    switch controller
        case 'jMPC';
            init0 =[m; s(:)];
            [m0, s0, predobj]= pgppred(dynmodel, m, s, u(:,i), 1);
            [ltvm.A, ltvm.B]= linearProbmodel(dynmodel, predobj.dMdm, ...
                predobj.dMds, predobj.dSdm, predobj.dSds);  
            ltvm.C = diag(ones(length(m)+length(s(:)), 1));  ltvm.D=0;   
            plantobj.A = ltvm.A;
            plantobj.B = ltvm.B;
            plantobj.C = ltvm.C;
            plantobj.D = ltvm.D;
            plantobj.Ts = 0.1;
    
            %     con.u = [-0.5,0.5, Inf; -1,1, Inf; -1,1, Inf];
            con.u = [];
            con.y = [];
            setp = [ref.x, repmat(0*s(:)', H, 1)];

            plantobj.con = con;
            plantobj.init0 = init0;
            controlobj.H = H;
            controlobj.ywt = ones(1, length(m)+length(s(:)));
            controlobj.uwt = [21000 27000 21000];
            controlobj.setpoints = setp;
    
            [optu, yp] = useJmpc(plantobj, controlobj);
        case 'mycontroller';
            [m0, s0, predobj]= pgppred(dynmodel, m, s, u(:,i), 1);
            [ltvm.A, ltvm.B]= linearProbmodel(dynmodel, predobj.dMdm, ...
                predobj.dMds, predobj.dSdm, predobj.dSds,m0,s0);  
            [dA, dB] = c2d(ltvm.A, ltvm.B, Ts);
            Q = diag(ones(dynmodel.sNum+dynmodel.sNum^2,1));
            R = diag([27000 21000 21000]);
            [tildeA, tildeB, tildeQ, tildeR] =...
                extendMatrix(H, dA, dB, Q, R);
            Bk = tildeB'*tildeQ*tildeB+ tildeR;
%             [HBk,exitflag] = chol(Bk);
%             if(exitflag)
%                 [Hes, jitter] = jit_chol(Bk);
%             else
%                 Hes = HBk;
%             end
            issym=@(x) all(all(tril(x)==triu(x).'));
            if(issym(Bk))
                Bk = Bk;
            else
                Bk = (Bk+Bk')/2;
            end
            Hes = Bk;
            setp = [ref.x, repmat(0*s(:)', H, 1)];
            
            x0 =[m; s(:)]-setp(1,:)';
            f = 2*(x0)'*tildeA'*tildeQ*tildeB;
            lbc = [umin(1).*ones(H,1)-ref.u(:,1),...
                umin(2).*ones(H,1)-ref.u(:,2),...
                umin(3).*ones(H,1)-ref.u(:,3)];
            ubc = [umax(1).*ones(H,1)-ref.u(:,1),...
                umax(2).*ones(H,1)-ref.u(:,2),...
                umax(3).*ones(H,1)-ref.u(:,3)];         
            lbc = reshape(lbc', H*dynmodel.iNum,1);
            ubc = reshape(ubc', H*dynmodel.iNum,1);
            
            optuu= quadprog(Hes,f,[],[],[],[],lbc,ubc,...
                deltaU-reshape(ref.u', H*dynmodel.iNum,1));
            optu= reshape(optuu', dynmodel.iNum,H)';
            optu = optu(1,:);
            yp=0;
    end
    
    %  linearised model around the operation point (x0,u0).
    u = [u, optu'+ref.u(1,:)'];
    [m, s]= pgppred(dynmodel, m, s, u(:,i+1), 1);
    Mmat(:,i) = m;
    Smat(:,:,i) =s;
    lMmat(:,i) = yp';

    fval(i) = sum((m([1 3 5])-ref.x(1,[1 3 5])').^2);    
    fprintf('LTVMPC:Iter %6.0d, StagError %12.8d\n', i, fval(i));
end
toc

res = Mmat';
lsyspred=lMmat';
reff = pathtargets(:,opt.out);

switch loadType
    case 'self0109';
        range1 = quad_ds.r1;
        range2 = quad_ds.r2;
        res = normalmatrix('reverse',res,range1(opt.out,:),range2);
        reff = normalmatrix('reverse',reff,range1(opt.out,:),range2);
        u = normalmatrix('reverse',u',range1(opt.in,:),range2);
        u=u';
    case 'Gauss01';
    otherwise;
end

gppred = res(:,[1 3 5]);
refpred = reff(:,[1 3 5]);
lsyspred = 0;

if nargout >1
    varargout{1} = refpred;
    varargout{2} = lsyspred;
    varargout{3} = u;
    varargout{4} = fval;
    varargout{5} = Smat;
end
end

function [A, B] = ...
    oldlinearProbmodel(dynmodel, dMdm, dMds, dSdm, dSds)
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

function [A, B] = ...
    linearProbmodel(dynmodel, dMdm, dMds, dSdm, dSds, m0, s0)
    D = dynmodel.iNum;
    E = dynmodel.sNum;
    
    dS0dsqrtS0 = 2*repmat(reshape(s0^(0.5), 1, numel(s0))', 1,numel(s0));
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

function [tileM, tileS] = jointDistribution(m, s, u)
    dimS = size(s, 1);
    dimU = size(u,1);
    tileM = [m; u];
    tileS = zeros(dimS+dimU,dimS+dimU);
    tileS(1:dimS, 1:dimS) = s;
end

function [optu, opty] = useJmpc(plantobj, controlobj)

    % jSS is based on using "jMPC" toolbox
    plant = jSS(plantobj.A, plantobj.B, plantobj.C, plantobj.D);
    plant = c2d(plant, plantobj.Ts);

    model = plant;
    plant.x0 =  plantobj.init0;
    model.x0 =  plantobj.init0;
    %     setp = zeros(1,42);

    %     %=====================
    %     % when  using specified QP solvers.
    %     %  'Matlab' --> all QP solvers
    %     %  'Mex'    -->  'wright' and 'mehrotra'
    %     % ------------------------------------------------------
    %      opts.QPSolver = 'mehrotra';
    %      MPC = jMPC(model, H, H, uwt,ywt,con,[],opts);
    %     %=====================
    
    MPC = jMPC(model, controlobj.H, controlobj.H,...
        controlobj.uwt, controlobj.ywt, plantobj.con);
    simopts = jSIM(MPC,plant, controlobj.H, controlobj.setpoints);
    simresult = sim(MPC,simopts,'Mex');
    optu =  simresult.plotvec.u(2,:);
    opty =  simresult.plotvec.yp(2,:);
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