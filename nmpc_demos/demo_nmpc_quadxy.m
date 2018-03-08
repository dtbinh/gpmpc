function res = demo_nmpc_quadxy

global ubar;
global eta;
ubar = 1;
eta = 0.2;

mpciterations = 200;
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
xr = path.xr;
yr = path.yr;
o = path.o;
trajectory = [xr', o, yr', o];
ref = trajectory(t,:);
Q = diag([27000 1 21000 1]);
R = diag([0.5 0.98]);
cost = (x-ref)*Q*(x-ref)' + u*R*u';
end

function cost = terminalcosts(t, x)
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

phi = pi/4;
theta = pi/4;
m = 0.74;
g = 9.81;
y = zeros(4, 1);
y(1) = x(2);
y(2) = u(1)/m;
y(3) = x(4);
y(4) = u(2)/m;
y = y';
end