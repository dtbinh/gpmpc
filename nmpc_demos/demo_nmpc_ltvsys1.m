function res = demo_nmpc_ltvsys1()

% External disturbances are not considered now
%    0.01*randn(1) may result in a big noise

clear;  clc; close all;

% [temp, path]=ode45('duffing',[0,40],[1.5;0]);

temp = pathGenerate(1, 201, 'Lorenz'); 
path = [temp.xr; temp.yr]';
assignin('base','path',path);

plot(path(:,1), path(:,2));
pause(3);

tmeasure = 0.0;
mpciterations = size(path, 1)-10;
N = 3;
T = 1;
xmeasure  = [0 0];
u0  = zeros(N, 2);

% path = TracjectorGenerator(1:210, 4);
% path = TracjectorGenerator(0:6*pi/210:6*pi, 2);

GNobj = GNgenerate(mpciterations+5, 2, 0, 1);
assignin('base','GNobj',GNobj);

[res.t, res.x, res.u] = mynmpc(@runningcosts, @terminalcosts, @constraints, ...
    @terminalconstraints, @linearconstraints, @system_dt, ...
    mpciterations, N, T, tmeasure, xmeasure, u0);

figure(1);
plot(res.x); hold on; plot(path,'k:');
xlim([0 mpciterations]);
% 
% obj.path = path;
% obj.x = res.x;
% obj.u = res.u;
% save lmpc4ltv_lorenz.mat obj;

ds = [res.x(1:end-1,:), res.u(1:end-1,:), res.x(2:end,:)];
% save lmpc4ltv_lorenz.txt ds -ascii
end

function cost = runningcosts(t, x, u)
    path = evalin('base','path');
    GNobj = evalin('base','GNobj');
    GNnoise = GNobj.UPlimit;
    Q = diag([1 1]);
    R = diag([1 1]);
    ref = path(t,:);
    cost = (x+0.01*GNnoise-ref)*Q*(x+0.01*GNnoise-ref)';
end

function cost = terminalcosts(t, x)
    cost =0;
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

GNobj = evalin('base','GNobj');
% a =  1+1*sin(2*pi*t/1500);
% b =  cos(2*pi*t/1500);

a = (1-exp(-t.^2))/(t+1);
b=1;
A = [1 0;0 1];
B = [a 1; b 0];
y = A*x'+B*u';
y = y' + 0.01*GNobj.GNnoise(t+1,1);

end
