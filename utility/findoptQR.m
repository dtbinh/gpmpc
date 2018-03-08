function [Q R] = findoptQR(system, t, xk, uk)

psoSwitchCG = 2; 
assignin('base','psoSwitchCG',psoSwitchCG);
% [Gbest,epoch,optParams] = psoLearn(@(params) costfunc(t, xk, uk, params));
[Gbest,epoch,optParams] = psoLearn(@(params) costfunc1(system, t, xk, params));

% Q = diag([optParams(1:4)]);
% R = diag([optParams(5:end)]);
Q = diag([optParams(1), 1, optParams(2), 1]);
R = diag([optParams(3), 1]);
assignin('base','updateQ',Q);
assignin('base','updateR',R);
end

function [cost] = costfunc(t, x, u, params)
Nc = size(params, 1);  % number of particles
Nx = size(x ,1);         
cost = zeros(1, Nc);
path = evalin('base','path');
xr = path.xr;
yr = path.yr;
o = path.o;
trajectory = [xr', o, yr', o];
sys = loadssfun_trans();

for i = 1:Nc 
    for j=1:Nx 
        ref = trajectory(t+j-1,:);
        Q = diag([params(i, 1:4)]);
        R = diag([params(i, 5:end)]);  
        cost(i) = cost(i)+ (x(j,:)-ref)*Q*(x(j,:)-ref)'+ u(j, :)*R*u(j, :)';
    end
    [K, P, E] = lqr(sys.Axy, sys.Bxy, Q, R);
    cost(i) = cost(i) + x(end,:)*P*x(end,:)';
end

% for i=1:Nc
%     Q = diag([params(i, 1), 1, params(i, 2), 1, params(i, 3), 1]);
%     R = diag([params(i, 4), 1]);  
% %     [K, P, E] = lqr(sys.Axy, sys.Bxy, Q, R);
% %    [K, P, E] = care(sys.Axy, sys.Bxy, Q, R);
% %    [K1, P1, E1] = dare(sys.Axy, sys.Bxy, Q, R);
%    [K2, P2, E2] = lqr(sys.Axy, sys.Bxy, Q, R);
%     for j=1:Nx
%         ref = trajectory(t+j,:);
%         cost(i) = cost(i)  + norm((sys.Axy-sys.Bxy*K2)*x(j,:)'-ref');
%     end
% end

end

function [cost] = costfunc1(system, t, x, params)

Nc = size(params, 1);
cost = zeros(1, Nc);
path = evalin('base','path');
xr = path.xr;
yr = path.yr;
o = path.o;
trajectory = [xr', o, yr', o];
sys = loadssfun_trans();

for i = 1:Nc 
    Q = diag([params(i, 1), 1, params(i, 2), 1]);
    R = diag([params(i, 3), 1]); 
    [K, P, E] = lqr(sys.Axy, sys.Bxy, Q, R);
    for j=1:3
        ref = trajectory(t+j-1,:);
        u0 = -K*x(j,:)';
        xx = system(t+j-1, x(j,:), u0', 1);
        cost(i) = cost(i) + norm(xx-ref);
    end
end
end

function [Gbest,epoch,optParams] = psoLearn(costfunc)

paraStruct = get_default_psoParams();
psoparams=struct2cell(paraStruct);
dNum = 3;
VarRange=repmat([1, 1000000], dNum-1, 1);
VarRange = [VarRange; [1, 100]];
[pso_out,tr,te] = ...
     pso_Trelea_vectorized(...
        costfunc, ...
        dNum,...
        4,...
        VarRange,0,...
        psoparams);
optParams = pso_out(1:end-1);
Gbest = te;
epoch= tr;

end