function res = nmpc_quad_rotate(mpciterations, anglepath)

global ubar;
global eta;
ubar = 1;
eta = 0.2;

assignin('base','anglepath',anglepath);
N = 3;
T = 1;
tmeasure = 0.0;
xmeasure  = [0 0 0 0 0 0];
u0  = zeros(N, 3);

[res.t, res.x, res.u, res.optobj] = mynmpc(@runningcosts, @terminalcosts, @constraints, ...
    @terminalconstraints, @linearconstraints, @system_dt, ...
    mpciterations, N, T, tmeasure, xmeasure, u0);

end

function cost = runningcosts(t, x, u)
    path = evalin('base','anglepath');
    GNobj = evalin('base','GNobj');
    % GNnoise = GNobj.UPlimit;
    GNnoise = [0, GNobj.UPlimit(1), ...
        0, GNobj.UPlimit(2),...
        0, GNobj.UPlimit(3)];
    xr = path.xr;
    yr = path.yr;
    zr = path.zr;
    o = path.o;
    trajectory = [xr, o, yr, o, zr, o];
    ref = trajectory(t,:);
    Q = diag([27000 1 21000 1 51000 1]);
    R = diag([1.3 1.5 1.1]);
    cost = (x+GNnoise-ref)*Q*(x+GNnoise-ref)' + u*R*u';
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

quad_params;
GNobj = evalin('base','GNobj');
y = zeros(1, 6);
% y(1) = x(2)+ GNobj.GNnoise(t+1, 1);
% y(2) = x(4)*x(6)*a1+x(4)*a2*OmegaR+b1*u(1)+ GNobj.GNnoise(t+1, 2);
% y(3) = x(4)+ GNobj.GNnoise(t+1, 3);
% y(4) = x(2)*x(6)*a3-x(2)*a4*OmegaR+b2*u(2)+ GNobj.GNnoise(t+1, 4);
% y(5) = x(6)+ GNobj.GNnoise(t+1, 5);
% y(6) = x(4)*x(2)*a5+b3*u(3)+ GNobj.GNnoise(t+1, 6);

y(1) = x(2);
y(2) = x(4)*x(6)*a1+x(4)*a2*OmegaR+b1*u(1)+ GNobj.GNnoise(t+1, 1);
y(3) = x(4);
y(4) = x(2)*x(6)*a3-x(2)*a4*OmegaR+b2*u(2)+ GNobj.GNnoise(t+1, 2);
y(5) = x(6);
y(6) = x(4)*x(2)*a5+b3*u(3)+ GNobj.GNnoise(t+1, 3);
end