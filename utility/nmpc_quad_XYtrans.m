function res = nmpc_quad_XYtrans(iters, path)

global ubar;
global eta;
ubar = 1;
eta = 0.2;

mpciterations = iters;
assignin('base','path',path);
N = 3;
T = 1;
tmeasure = 0.0;
xmeasure  = [0 0 0 0];
u0  = zeros(N, 2);

[res.t, res.x, res.u] = mynmpc(@runningcosts, @terminalcosts, @constraints, ...
    @terminalconstraints, @linearconstraints, @system_dt, ...
    mpciterations, N, T, tmeasure, xmeasure, u0);
end

function cost = runningcosts(t, x, u)
path = evalin('base','path');
GNobj = evalin('base','GNobj');
GNnoise = GNobj.UPlimit(1:4);
xr = path.xr;
yr = path.yr;
o = path.o;
trajectory = [xr', o, yr', o];
ref = trajectory(t,:);
sys = loadssfun_trans();

if ~isempty(evalin('base','updateQ'))
    updateQ = evalin('base','updateQ');
    updateR = evalin('base','updateR');
    cost = (x+GNnoise-ref)*updateQ*(x+GNnoise-ref)'...
        + u*updateR*u';
else
    cost = (x+GNnoise-ref)*sys.Qxy*(x+GNnoise-ref)' + u*sys.Rxy*u';
end

end

function cost = terminalcosts(t, x)
% path = evalin('base','path');
% GNobj = evalin('base','GNobj');
% GNnoise = GNobj.UPlimit(1:4);
% xr = path.xr;
% yr = path.yr;
% o = path.o;
% trajectory = [xr', o, yr', o];
% ref = trajectory(t,:);
% sys = loadssfun_trans();
% if ~isempty(evalin('base','updateQ'))
%     updateQ = evalin('base','updateQ');
%     updateR = evalin('base','updateR');
%     [K, P, E] = lqr(sys.Axy, sys.Bxy, updateQ, updateR);
%     cost =  (x+GNnoise-ref)*P*(x+GNnoise-ref)';
% else
%     [K, P, E] = lqr(sys.Axy, sys.Bxy, sys.Qxy, sys.Rxy);
%     cost =  (x+GNnoise-ref)*P*(x+GNnoise-ref)';
% end
cost = 0.0;
end
function [c,ceq] = constraints(t, x, u)
    global eta;
    c(1) =   Inf;
    c(2) =  -Inf;
    c(3) =  Inf;
    c(4) =  -Inf;
    c(5) =  Inf;
    c(6) =  -Inf;
    c(7) =  Inf;
    c(8) =  -Inf;
    ceq  = [];
end
function [c,ceq] = terminalconstraints(t, x)
    global eta;
    c(1) =   Inf;
    c(2) =  -Inf;
    c(3) =   Inf;
    c(4) =  -Inf;
    c(5) =   Inf;
    c(6) =  -Inf;
    c(7) =   Inf;
    c(8) =  -Inf;
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
% GNobj = GNgenerate(1, 4, 0, 0.1);
sys = loadssfun_trans();
y = sys.Axy*x'+sys.Bxy*u';
y = y' + GNobj.GNnoise(t+1, 1:4);

end