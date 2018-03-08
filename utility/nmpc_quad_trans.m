function res = nmpc_quad_trans(iters, path, varargin)
%----------------------------------------------------------------
%   varargin{1}  'none',   no multimedia generated
%----------------------------------------------------------------

mpciterations = iters;
assignin('base','path',path);
N = 3;
T = 1;
tmeasure = 0.0;
if nargin>2 && ~isempty(varargin{1})
    xmeasure = varargin{1};
else
    xmeasure  = [0 0 0 0 0 0];
end
u0  = zeros(N, 3);

[res.t, res.x, res.u, res.optobj] = mynmpc(@runningcosts, @terminalcosts, @constraints, ...
    @terminalconstraints, @linearconstraints, @system_dt, ...
    mpciterations, N, T, tmeasure, xmeasure, u0);
end
function cost = runningcosts(t, x, u)
path = evalin('base','path');
GNobj = evalin('base','GNobj');
% GNnoise = GNobj.UPlimit;
GNnoise = [0, GNobj.UPlimit(1), ...
    0, GNobj.UPlimit(2),...
    0, GNobj.UPlimit(3)];
xr = path.xr;
yr = path.yr;
zr = path.zr;
o = path.o;
trajectory = [xr', o, yr', o, zr', o];
ref = trajectory(t,:);

refu = zeros(1, 3);
refu(3) = 0.74*(ref(3)+9.81);
refu(1) = 0.74*ref(1)/refu(3);
refu(2) = 0.74*ref(2)/refu(3);

%  These parameters are good for system without noises
Q = diag([210000 1 270000 1 210000 1]);
R = diag([1 1 1]);
% cost = (x+GNnoise-ref)*Q*(x+GNnoise-ref)' ...
%     + (u-refu)*R*(u-refu)';

cost = (x-ref)*Q*(x-ref)' ...
    + u*R*u';
end
function cost = terminalcosts(t, x)
    P = diag([0 0 0 0 0 0]);
    cost = x*P*x';
end
function [c,ceq] = constraints(t, x, u)
    global eta;
    c(1) =   Inf;
    c(2) =  -Inf;
    c(3) =  Inf;
    c(4) =  -Inf;
    c(5) =  Inf;
    c(6) =  -Inf;
    c(7) =   Inf;
    c(8) =   Inf;
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
%     windgust = windspeed(t,50,120,'step',0,2);
    windgust = 0*windspeed(t,50,120,'sine',0,2);
    GNobj = evalin('base','GNobj');
    phi = pi/4;  theta = pi/4;
    m = 0.74;    g = 9.81;
    y = zeros(6, 1);
%     y(1) = x(2) + GNobj.GNnoise(t+1, 1);
%     y(2) = u(2)*u(1)/m + GNobj.GNnoise(t+1, 2);
%     y(3) = x(4)+ GNobj.GNnoise(t+1, 3);
%     y(4) = u(3)*u(1)/m+ GNobj.GNnoise(t+1, 4);
%     y(5) = x(6)+ GNobj.GNnoise(t+1, 5);
%     y(6) = g-(cos(phi)*cos(theta))*u(1)/m+ GNobj.GNnoise(t+1, 6);
%     y(6) = g-(cos(t/4)*cos(t/2))*u(1)/m+ GNobj.GNnoise(t+1, 6);

    y(1) = x(2);
    y(2) = u(2)*u(1)/m + GNobj.GNnoise(t+1, 1)+windgust;
    y(3) = x(4);
    y(4) = u(3)*u(1)/m+ GNobj.GNnoise(t+1, 2);
    y(5) = x(6);
    y(6) = g-(cos(phi)*cos(theta))*u(1)/m+ GNobj.GNnoise(t+1, 3);
    y = y';
end

