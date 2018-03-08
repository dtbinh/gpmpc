clc; clear; close all;
% rmpath(genpath('C:\Users\gcao\Desktop\mpt'));

% Lorenz   Elliptical   Helix   crown   spraypool  randpath
pathtype = 'Lorenz';    
mpciterations = 199;   %  180    199    89

% fminsearch  fminlbfgs  fmincon   minimize
optimizer = 'fmincon';
[t_gppred, t_syspred, t_refpred, t_u, M1, S1, fval1] = ...
    pgp_smpc_quad(pathtype, mpciterations, 'TRAN', optimizer);

[r_gppred, r_syspred, r_refpred, r_u, M2, S2, fval2] = ...
    pgp_smpc_quad(pathtype, mpciterations, 'ROTATE', optimizer);

paraset = getRouteParams(pathtype);
plotvec_t = 1:mpciterations-10;
plotvec_y = t_gppred;
plotvec_ref = t_refpred(plotvec_t,:);
plotvec_angle = r_gppred;

Tobj.gppred = t_gppred;
Tobj.refpred = t_refpred;
Tobj.u = t_u;
Tobj.S1 = S1;

Robj.gppred = r_gppred;
Robj.refpred = r_refpred;
Robj.u = r_u;
Robj.S1 = S2;

quadTranRotaplot3(plotvec_t, plotvec_y, plotvec_ref, ...
            plotvec_angle, 'none', [], ...
            paraset.flyerpara,paraset.viewangle);

% figure(2);
% drawRealPrediction(plotvec_ref,plotvec_y,0);
% legend('Trajectory', 'GP Mean');
% 
% figure(3);
% plotvec_ref = r_refpred(plotvec_t,:);
% drawRealPrediction(plotvec_ref,plotvec_angle,0);
% legend('Trajectory', 'GP Mean');
% 
% figure(4);
% cu = t_u';
% subplot(2,1,1); plot(cu(:,1), 'k-');  legend('u_1');
% subplot(2,1,2); plot(cu(:,2), 'b-'); hold on;
% plot(cu(:,3), 'k-'); legend('u_x', 'u_y');
% 
% figure(5);
% cu = r_u';
% subplot(2,1,1); plot(cu(:,1), 'k-'); hold on;
% plot(cu(:,2), 'b-'); legend('u_2', 'u_3');
% subplot(2,1,2); plot(cu(:,3), 'k-');legend( 'u_4');
% 
% figure(6);
% subplot(1,2,1); plot(fval1); 
% ylabel('Cost Values'); xlabel('Steps');
% subplot(1,2,2); plot(fval2); xlabel('Steps');


% obj.gppred = t_gppred;
% obj.refpred = t_refpred;
% obj.u = t_u';
% for i=1:170
%     Smat(:,:,i) =S1{i};
% end
% obj.Smat = Smat;
% save Tobj_fminsearch_lorenz.mat obj;
% 
% clear obj;
% obj.gppred = r_gppred;
% obj.refpred = r_refpred;
% obj.u = r_u';
% for i=1:170
%     Smat(:,:,i) =S2{i};
% end
% obj.Smat = Smat;
% save Robj_fminsearch_lorenz.mat obj;

