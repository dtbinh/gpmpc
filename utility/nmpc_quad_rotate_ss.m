function res = nmpc_quad_rotate_ss(mpciterations, anglepath)

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

[res.t, res.x, res.u, res.optobj] = nmpc(@runningcosts, @terminalcosts, @constraints, ...
    @terminalconstraints, @linearconstraints, @system_dt, ...
    mpciterations, N, T, tmeasure, xmeasure, u0);

end

function cost = runningcosts(t, x, u)
    path = evalin('base','anglepath');
    xr = path.xr;
    yr = path.yr;
    zr = path.zr;
    o = path.o;
    trajectory = [xr, o, yr, o, zr, o];
    ref = trajectory(t,:);
    sys = loadssfun_rotate();
    cost = (x-ref)*sys.Q*(x-ref)' + u*sys.R*u';
end

function cost = terminalcosts(t, x)
    sys = loadssfun_rotate();
    [K, P, E] = lqr(sys.A, sys.B, sys.Q, sys.R);
    cost =  x*P*x';
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

% quad_params;
% y = zeros(1, 6);
% y(1) = x(2);
% y(2) = x(4)*x(6)*a1+x(4)*a2*OmegaR+b1*u(1);
% y(3) = x(4);
% y(4) = x(2)*x(6)*a3-x(2)*a4*OmegaR+b2*u(2);
% y(5) = x(6);
% y(6) = x(4)*x(2)*a5+b3*u(3);
scale = 0.1;
noise=scale*randn(6,1);
% noise(1) = 0;noise(3) = 0;noise(5) = 0;
sys = loadssfun_rotate();
y = sys.A*x'+sys.B*u' + noise;
y = y';
end