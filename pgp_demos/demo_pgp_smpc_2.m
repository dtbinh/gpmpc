clc; clear; close all;
pathtype = 'Lorenz';    % Lorenz   Elliptical   Helix   crown   spraypool
mpciterations = 199;   

filename1 = strcat('ds_nmpc_', pathtype, '_trans');
loadType = 'self0109';  % original self0109

tempRandSequences = randperm(mpciterations);
opt.dsSampleRange = sort(tempRandSequences(1:mpciterations));
opt.dsTestRange = 2:mpciterations;

% %  trans_z control
% opt.state = [5:7];
% opt.in = 7;
% opt.out = 14:15;
% 
% quad_ds = loadSampleData(filename1,loadType,opt);
% % numtem = numel(quad_ds);
% 
% dynmodel.loadType = loadType;
% dynmodel.ds = quad_ds;
% dynmodel.inputs = quad_ds.proData(opt.dsSampleRange,opt.state);
% dynmodel.targets = quad_ds.proData(opt.dsSampleRange,opt.out);
% dynmodel.train = @train;
% dynmodel. iNum =1;
% dynmodel. sNum =2;
% dynmodel. oNum =2; 
% 
% mpcobj.N = mpciterations-10;
% mpcobj.H = 3;
% 
% [tz_gppred, tz_refpred, tz_u, Mz, Sz, fvalz] = ...
%     pgp_smpc(dynmodel, mpcobj, 'TRAN_Z');

%  trans_xy control
opt.state = [1:4, 8:10];
opt.in = 8:10;
opt.out = 10:13;

quad_ds = loadSampleData(filename1,loadType,opt);
% numtem = numel(quad_ds);

dynmodel.loadType = loadType;
dynmodel.ds = quad_ds;
dynmodel.inputs = quad_ds.proData(opt.dsSampleRange,opt.state);
dynmodel.targets = quad_ds.proData(opt.dsSampleRange,opt.out);
dynmodel.train = @train;
dynmodel. iNum =3;
dynmodel. sNum =4;
dynmodel. oNum =4; 

mpcobj.N = mpciterations-10;
mpcobj.H = 3;

[txy_gppred, txy_refpred, txy_u, Mxy, Sxy, fvalxy] = ...
    pgp_smpc(dynmodel, mpcobj, 'TRAN_XY');

t_gppred = [txy_gppred, tz_gppred];
t_refpred = [txy_refpred, tz_refpred];
t_u = [tz_u; txy_u];
%---------------------------------------------------------
filename2 = strcat('ds_nmpc_', pathtype, '_rotate');

tempRandSequences = randperm(mpciterations);
opt.dsSampleRange = sort(tempRandSequences(1:mpciterations));
opt.dsTestRange = 2:mpciterations;

opt.state = 1:9;
opt.in = 7:9;
opt.out = 10:15;
quad_ds = loadSampleData(filename1,loadType,opt);

dynmodel.ds = quad_ds;
dynmodel.inputs = quad_ds.proData(opt.dsSampleRange,opt.state);
dynmodel.targets = quad_ds.proData(opt.dsSampleRange,opt.out);
dynmodel.train = @train;
dynmodel. iNum =3;
dynmodel. sNum =6;
dynmodel. oNum =6; 

mpcobj.N = mpciterations-10;
mpcobj.H = 1;

[r_gppred, r_refpred, r_u, M2, S2, fval2] = ...
    pgp_smpc(dynmodel, mpcobj, 'ROTATE');

paraset = getRouteParams(pathtype);
plotvec_t = 1:mpciterations-10;
plotvec_y = t_gppred(:,[1 3 5]);
plotvec_ref = t_refpred(plotvec_t,[1 3 5]);
plotvec_angle = r_gppred(:, [1 3 5]);

quadTranRotaplot3(plotvec_t, plotvec_y, plotvec_ref, ...
            plotvec_angle, 'none', [], ...
            paraset.flyerpara,paraset.viewangle);

figure(2);
drawRealPrediction(plotvec_ref,plotvec_y,0);
legend('Trajectory', 'GP Mean');

figure(3);
cu = t_u';
subplot(2,1,1); plot(cu(:,1), 'k-');  legend('u_1');
subplot(2,1,2); plot(cu(:,2), 'b-'); hold on;
plot(cu(:,3), 'k-'); legend('u_x', 'u_y');

figure(4);
cu = r_u';
subplot(2,1,1); plot(cu(:,1), 'k-'); hold on;
plot(cu(:,2), 'b-'); legend('u_2', 'u_3');
subplot(2,1,2); plot(cu(:,3), 'k-');legend( 'u_4');

figure(5);
subplot(1,3,1); plot(fvalxy); 
ylabel('Cost Values'); xlabel('Steps');
subplot(1,3,2); plot(fvalz); xlabel('Steps');
subplot(1,3,3); plot(fval2); xlabel('Steps');


