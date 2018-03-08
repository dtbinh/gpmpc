%%
% Quadrotor fully control using nonlear mpc strategy 
%  without consideration of external disturbances
%       (x-ref)'*Q*(x-ref) + u'*R*u + terminal cost;
%%

clear all; close all; clc;

%  Lorenz   Elliptical   Helix   crown   spraypool   randpath
pathtype = 'Lorenz';  

switch pathtype
    case 'switchpath';
        mpciterations = 399;
        tempt = 0:mpciterations*2;
        path = BustrophedonPath_Generate(numel(tempt));
        path.o = zeros(numel(tempt), 1);
        axr = [1 6 1 5 -0.5 1.5];
    case {'SingleRing', 'Hover', 'Elliptical', 'Helix', 'Lorenz'};
        mpciterations = 200;
        path = pathGenerate(0.5, mpciterations+1, pathtype);
        axr = [];
    case 'crown';
        mpciterations = 200;
        t=-pi:2*pi/(mpciterations+10):pi;
        path.xr = cos(t);
        path.yr = sin(t);
        path.zr = sin(5*t);
        path.o = zeros(numel(t),1);
        path.t = t;      
        axr = [-2 2 -2 2 -2 2];
    case'spraypool';
        mpciterations = 128;
        t = mpciterations + 10;
        th = (1:t)/t*2*pi;
        path.xr = cos(th);
        path.yr = sin(th);
        path.zr = abs(fft(ones(10,1),t))';
        path.o = zeros(t,1);
        path.t = th;
        axr = [-2 2 -2 2 0 20];
    case'randpath';
        mpciterations = 90;
        t=1:99;
        temp = load('randpath.txt');
        path.xr = temp(t,1)';
        path.yr = temp(t,2)';
        path.zr = temp(t,3)';
        path.o = zeros(max(t),1);
        path.t =t;
        axr = [-1 1 -1 1 0 1.5];
    otherwise;
        error('');
end
% txtfilename = strcat('ds_nmpc_', pathtype);
% gpobj = loadgptrainobj(txtfilename);
% assignin('base','GPobj',gpobj);

GNobj = GNgenerate(mpciterations+5, 3, 0, 1);
assignin('base','GNobj',GNobj);
res = nmpc_quad_trans(mpciterations+5, path);
% nmpc_trans_plot(res, mpciterations, path);
% plot(res.optobj.iteration, res.optobj.objvalue)

GNobj = GNgenerate(mpciterations+5, 3, 0, 1);
assignin('base','GNobj',GNobj);
ref_phi = -asin(res.u(:,3));
ref_theta = asin(1./cos(angle(ref_phi)).*res.u(:,2));
% ref_psi = zeros(numel(ref_phi), 1)+0.53;
ref_psi = zeros(numel(ref_phi), 1)+0;
o = zeros(numel(ref_phi),1);
anglepath.xr = ref_phi;
anglepath.yr = ref_theta;

anglepath.zr = ref_psi;
anglepath.o = o;

angleres = nmpc_quad_rotate(mpciterations, anglepath);
% nmpc_rotate_plot(angleres, mpciterations, anglepath);
%plot(angleres.optobj.iteration, angleres.optobj.objvalue)

paraset = getRouteParams(pathtype);
plotvec_t = 1:mpciterations;
plotvec_y = [res.x(plotvec_t,1), ...
    res.x(plotvec_t,3), res.x(plotvec_t,5)];
plotvec_ref = [path.xr(plotvec_t)', ...
    path.yr(plotvec_t)', path.zr(plotvec_t)'];
plotvec_angle = [angleres.x(plotvec_t,1), ...
    angleres.x(plotvec_t,3), angleres.x(plotvec_t,5)];

quadTranRotaplot3(plotvec_t, plotvec_y, plotvec_ref, ...
            plotvec_angle, 'none', axr, ...
            paraset.flyerpara,paraset.viewangle);

% txt_ds_trans = [res.x(1:end-1,:), res.u(1:end-1,:), ...
%     res.x(2:end, :)];
% txt_ds_rotate = [angleres.x(1:end-1,:),...
%     angleres.u(1:end-1,:), angleres.x(2:end, :)];
% tempname1 = strcat('ds_nmpc_', pathtype, '_trans.txt');
% tempname2 = strcat('ds_nmpc_', pathtype, '_rotate.txt');
% save(tempname1, 'txt_ds_trans', '-ascii');
% save(tempname2, 'txt_ds_rotate', '-ascii');
