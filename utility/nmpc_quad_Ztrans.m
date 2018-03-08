function res = nmpc_quad_Ztrans(iters, path)

global ubar;
global eta;
ubar = 1;
eta = 0.2;

mpciterations = iters;
assignin('base','path',path);
N = 3;
T = 1;
tmeasure = 0.0;
xmeasure  = [0 0];
u0  = zeros(N, 1);

[res.t, res.x, res.u] = mynmpc(@runningcosts, @terminalcosts, @constraints, ...
    @terminalconstraints, @linearconstraints, @system_dt, ...
    mpciterations, N, T, tmeasure, xmeasure, u0);

end

function cost = runningcosts(t, x, u)
path = evalin('base','path');
GNobj = evalin('base','GNobj');
GNnoise = GNobj.UPlimit(5:6);
zr = path.zr;
o = path.o;
trajectory = [zr', o];
ref = trajectory(t,:);
sys = loadssfun_trans();
cost = (x+GNnoise-ref)*sys.Qz*(x+GNnoise-ref)' + u*sys.Rz*u';
end

function cost = terminalcosts(t, x)
%     sys = loadssfun_trans();
%     path = evalin('base','path');
%     GNobj = evalin('base','GNobj');
%     GNnoise = GNobj.UPlimit(5:6);
%     zr = path.zr;
%     o = path.o;
%     trajectory = [zr', o];
%     ref = trajectory(t,:);
%     [K, P, E] = lqr(sys.Az, sys.Bz, sys.Qz, sys.Rz);
%     cost =  (x+GNnoise-ref)*P*(x+GNnoise-ref)';
cost = 0.0;
end
function [c,ceq] = constraints(t, x, u)
    global eta;
    c(1) =   Inf;
    c(2) =  Inf;
    c(3) =  Inf;
    c(4) =  Inf;
    c(5) =  Inf;
    c(6) =  Inf;
    ceq  = [];
end
function [c,ceq] = terminalconstraints(t, x)
    global eta;
    c(1) =   Inf;
    c(2) =  Inf;
    c(3) =   Inf;
    c(4) =  Inf;
    c(5) =   Inf;
    c(6) =  Inf;
    ceq  = [];
end
function [A, b, Aeq, beq, lb, ub] = linearconstraints(t, x, u)
    global ubar
    A   = [];
    b   = [];
    Aeq = [];
    beq = [];
    lb  = -Inf;
    ub  =  Inf;
end
function y = system_dt(t, x, u, T)

GNobj = evalin('base','GNobj');
sys = loadssfun_trans();
y = sys.Az*x'+sys.Bz*u' ;
y = y' + GNobj.GNnoise(t+1, 5:6);

end

