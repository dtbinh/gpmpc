function nmpc_rotate_plot(res, mpciterations, path)

clf; 
figure(1);
subplot(3,1,1); plot(res.t, res.x(:,1), 'b-'); hold on;
plot(res.t, path.xr(1:mpciterations), 'k--');
subplot(3,1,2); plot(res.t, res.x(:,3), 'b-'); hold on;
plot(res.t, path.yr(1:mpciterations), 'k--');
subplot(3,1,3); plot(res.t, res.x(:,5), 'b-'); hold on;
plot(res.t, path.zr(1:mpciterations), 'k--');
legend('Controlled', 'Reference');

figure(2);
subplot(3,1,1); plot(res.t, res.u(:,1),'k-'); 
legend('u_1'); hold on;
subplot(3,1,2); plot(res.t, res.u(:,2),'b-');
legend('u_2'); hold on;
subplot(3,1,3); plot(res.t, res.u(:,3),'r-');
legend('u_3');

return