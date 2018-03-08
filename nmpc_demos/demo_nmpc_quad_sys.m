%%
%  Nonlinear MPC for Quadrotor Control Demo
%  Translational system is divided into xy and z;
%%
clear all; close all; clc;

% 'SingleRing'   'Hover'   'Elliptical'   'Helix'  'Lorenz'
path = pathGenerate(0.5, 300, 'Elliptical');
assignin('base','path',path);

resxy = demo_nmpc_quadxy;
resz = demo_nmpc_quadz;
resvec.t = resxy.t;
resvec.x = [resxy.x(:,1), resxy.x(:,3), resz.x(:,1)];
resvec.u = [resxy.u, resz.u];

figure(1);
subplot(3,1,1); plot(resvec.t, resvec.x(:,1), 'b-'); hold on;
plot(resvec.t, path.xr(1:200), 'k--');
subplot(3,1,2); plot(resvec.t, resvec.x(:,2), 'b-'); hold on;
plot(resvec.t, path.yr(1:200), 'k--');
subplot(3,1,3); plot(resvec.t, resvec.x(:,3), 'b-'); hold on;
plot(resvec.t, path.zr(1:200), 'k--');

figure(2);
plot(resvec.t, resvec.u(:,1),'k-'); hold on;
plot(resvec.t, resvec.u(:,2),'b-'); hold on;
plot(resvec.t, resvec.u(:,3),'r-'); 

figure(3);
plot3(path.xr(1:200), path.yr(1:200), path.zr(1:200), 'r--'); hold on;
plot3(resvec.x(:,1), resvec.x(:,2), resvec.x(:,3), 'k-');
