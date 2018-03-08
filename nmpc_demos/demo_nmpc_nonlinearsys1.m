function res = demo_nmpc_nonlinearsys1()

% External disturbances are not considered now
%    0.01*randn(1) may result in a big noise

clear all;  close all; clc;
global lbset ubset;
tmeasure = 0.0;

N = 3;
T = 1;
xmeasure  = [0 0 0 0];
u0  = zeros(N, 2);

pathcase = 'lorenz';   %   lorenz,   sincos,  duffing
switch pathcase
    case 'duffing';
        mpciterations = 200;
        [temp, path]=ode45('duffing',[0,40],[1.5;0]);
    case 'lorenz';
        mpciterations = 200;
        temp = pathGenerate(1, mpciterations+10, 'Lorenz'); 
        path = [temp.xr; temp.yr]';
        lbset = [-4 -5];
        ubset = [4 5];
    case 'sincos';
        mpciterations = 175;
        path = TracjectorGenerator(0:6*pi/210:6*pi, 2);
    case 'step';
        mpciterations = 200;
        path = TracjectorGenerator(1:mpciterations+10, 4);
end

assignin('base','path',path);

GNobj = GNgenerate(mpciterations+5, 2, 0, 1);
assignin('base','GNobj',GNobj);

[res.t, res.x, res.u, res.optobj] = mynmpc(@runningcosts, @terminalcosts, @constraints, ...
    @terminalconstraints, @linearconstraints, @system_dt, ...
    mpciterations, N, T, tmeasure, xmeasure, u0);

subplot(2,1,1);
plot(res.x(:,[1 3])); hold on; plot(path,'k:'); xlim([0 mpciterations]);
subplot(2,1,2);
plot(res.u);

for i=1:numel(res.t)-1
    syssim(i,:) = system_dt(res.t(i+1), res.x(i,:), res.u(i,:));
end

res.syssim = syssim;

for i=1:numel(res.t)
    syssim0(i,:) = system_dt(res.t(i), res.x(i,:), res.u(i,:));
end

figure(1);
for i=1:4
    subplot(4,1,i);
    plot(res.x(:,i), 'k:'); hold on; plot(res.syssim(:,i),'b-');
    hold on; plot(syssim0(:,i), 'r-');
    xlim([1 numel(res.t)-1]);
end


% ds=[res.x(1:end-1,:), res.u(1:end-1,:), res.x(2:end,:)];
% tempname = strcat('ds_nmpc_nltv_', pathcase, '.txt');
% save(tempname, 'ds', '-ascii');
% 
% obj.path = path;
% obj.x = res.x;
% obj.u = res.u;
% tempname = strcat('nmpc_nltv_', pathcase, '.mat');
% save(tempname, 'obj');
end

function cost = runningcosts(t, x, u)
    path = evalin('base','path');
    GNobj = evalin('base','GNobj');
    GNnoise = GNobj.UPlimit;
    Q = diag([27000  21000]);
    R = diag([1 1]);
    ref = path(t,:);
    cost = (x([1 3])-ref)*Q*(x([1 3])-ref)'+...
        u*R*u';
end

function cost = terminalcosts(t, x)
    cost =0;
end

function [c,ceq] = constraints(t, x, u)
    global lbset ubset;
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
    global lbset ubset;
    c(1) =   Inf;
    c(2) =  -Inf;
    c(3) =   Inf;
    c(4) =  -Inf;
    c(5) =   Inf;
    c(6) =  -Inf;
    ceq  = [];
end

function [A, b, Aeq, beq, lb, ub] = linearconstraints(t, x, u)
    global lbset ubset;
    A   = [];
    b   = [];
    Aeq = [];
    beq = [];
    lb  = lbset;
    ub  =  ubset;
end

function y = system_dt(t, x, u, T)

GNobj = evalin('base','GNobj');

% a = 1 + 0.1*sin(2*pi*t/1500);
% b = 1 + 0.1*cos(2*pi*t/1500);

a =  10+0.5*sin(t);
b=  10./(1+exp(-0.05*t))-0.5;

y(1) = x(1)^2/(1+x(1)^2)+0.3*x(2)+0*0.01*GNobj.GNnoise(t+1, 1);
y(2) = x(1)^2/(1+x(2)^2+x(3)^2+x(4)^2) + b*u(1);
y(3) = x(3)^2/(1+x(3)^2)+0.2*x(4)+0*0.01*GNobj.GNnoise(t+1, 2);
y(4) = x(3)^2/(1+x(1)^2+x(2)^2+x(4)^2) + a*u(2);

end
