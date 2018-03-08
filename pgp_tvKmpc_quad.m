function [gppred, refpred, varargout] = ...
    pgp_tvKmpc_quad(pathtype, mpciterations, systype)
%%
filename1 = strcat('ds_nmpc_', pathtype, '_trans');
filename2 = strcat('ds_nmpc_', pathtype, '_rotate');

switch systype
    case 'TRAN';
        filename = filename1;
        loadType = 'original';  % original self0109
        H = 5; 
    case 'ROTATE';
        filename = filename2;
        loadType = 'self0109';
        H = 1;
    otherwise;
end

pilcoSettings;

tic
% 1. training dynamic GP model
[dynmodel nlml] = dynmodel.train(dynmodel, plant, [300 500]);
process_mes = strcat('Successfully Training', 32,  systype, ' GP Model\n');
fprintf(process_mes);
toc

% 2. time-varying linear policy K learning  (u=Kx);
pathtargets = quad_ds.proData;
N = 30;
mu0Sim = zeros(6, N+1); mu0Sim(:, 1) = [0 0 0 0 0 0]';
S0Sim = zeros(6,6,N+1);  S0Sim(:,:,1) = diag([0.1 0.1 0.1 0.1 0.1 0.1].^2);

fprintf('SMPC Starts...\n');
for i =1:N
    opt.fh = 1;
    traject = pathtargets(1:H, plant.dyno);
    
    %%
%     [policy.p fval] = minimize(policy.p, 'value', optparam, mu0Sim(:, i),...
%         S0Sim(:,:,i), dynmodel, policy, plant, cost, H, traject);

    %%
    %  nvalue function minimization using gradien-based methods
    %  here policy.p structure is reformed as a numeric vector.
    %active-set   sqp  interior-point  trust-region-reflective
    % %----------------------------------------------------------------
    options = optimoptions('fmincon','Algorithm','trust-region-reflective', 'GradObj', 'on'); 
        p0 = policy.p;
        u0 = unwrap(p0);
        Ulimit = zeros(1, numel(u0))+Inf;
        [optp, fval] = fmincon(@(p) nvalue(p, mu0Sim(:, i), S0Sim(:,:,i), dynmodel, policy,...
            plant, cost, H, traject), u0, [], [], [], [], -Ulimit, Ulimit, [], options);
%%
    toc
    [tempmu0 tempS0] = ...
        pred(policy, plant, dynmodel, mu0Sim(:, i), S0Sim(:,:,i), H);
    mu0Sim(:,i+1) = tempmu0(:,end);
    S0Sim(:,:,i+1) = tempS0(:,:,end);
    Policy{i}.p.w = policy.p.w;
    Policy{i}.p.b = policy.p.b;
    [um(:,i) us(:,:,i)] = conlin(policy, mu0Sim(:,i), S0Sim(:,:,i));
    fprintf('SMPC:Iter %6.0d , Fval %12.8d\n', i, fval(end));
end
toc
reff = pathtargets(:,opt.out);

switch loadType
    case 'self0109';
        range1 = quad_ds.r1;
        range2 = quad_ds.r2;
        mu0Sim = normalmatrix('reverse',mu0Sim',range1(opt.out, :),range2);
        reff = normalmatrix('reverse',reff,range1(opt.out,:),range2);
        um = normalmatrix('reverse',um,range1(opt.in,:),range2);
    case 'Gauss01';
    otherwise;
end

gppred = [res(:,1), res(:,3), res(:,5)];
refpred = [reff(:,1), reff(:,3), reff(:,5)];

if nargout >2
    varargout{1} = um;
    varargout{2} = S0Sim;
    varargout{3} = us;
    varargout{4} = fval;
    varargout{5} = Policy;
end
end








