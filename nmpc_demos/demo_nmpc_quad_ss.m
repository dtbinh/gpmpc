%%
% Quadrotor fully control using nonlear mpc strategy 
%  with consideration of external disturbances
%       (x-ref)'*Q*(x-ref) + u'*R*u + terminal cost;
%   Terminal cost = x'*P*x;    
%    s.t.       P-(A+B*K)'*P*(A+BK)-Q-K'*R*K >0
%%

clear all; close all; clc;
pathtype = 'Elliptical';
Q =[]; R=[];
assignin('base','updateQ',Q);
assignin('base','updateR',R);
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
    otherwise;
end

GNobj = GNgenerate(mpciterations+10, 6, 0, 0.1);
assignin('base','GNobj',GNobj);

% optModel = load('optModel_elliptical.mat');
resxy = nmpc_quad_XYtrans(mpciterations+5, path);
resz = nmpc_quad_Ztrans(mpciterations+5, path);
res.t = resxy.t;
res.x = [resxy.x, resz.x];
res.u = [ resz.u, resxy.u];

%nmpc_trans_plot(res, mpciterations, path);
%plot(res.optobj.iteration, res.optobj.objvalue)

% ref_phi = -asin(res.u(:,3));
% ref_theta = asin(1./cos(angle(ref_phi)).*res.u(:,2));
% % ref_psi = zeros(numel(ref_phi), 1)+0.53;
% ref_psi = zeros(numel(ref_phi), 1)+0;
% o = zeros(numel(ref_phi),1);
% anglepath.xr = ref_phi;
% anglepath.yr = ref_theta;
% anglepath.zr = ref_psi;
% anglepath.o = o;
% 
% angleres = nmpc_quad_rotate(mpciterations, anglepath);
% % nmpc_rotate_plot(angleres, mpciterations, anglepath);
% %plot(angleres.optobj.iteration, angleres.optobj.objvalue)
% 
% paraset = getRouteParams(pathtype);
% plotvec_t = 1:mpciterations;
% plotvec_y = [res.x(plotvec_t,1), ...
%     res.x(plotvec_t,3), res.x(plotvec_t,5)];
% plotvec_ref = [path.xr(plotvec_t)', ...
%     path.yr(plotvec_t)', path.zr(plotvec_t)'];
% plotvec_angle = [angleres.x(plotvec_t,1), ...
%     angleres.x(plotvec_t,3), angleres.x(plotvec_t,5)];
% 
% quadTranRotaplot3(plotvec_t, plotvec_y, plotvec_ref, ...
%             plotvec_angle, 'none', axr, ...
%             paraset.flyerpara,paraset.viewangle);


