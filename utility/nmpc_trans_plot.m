function nmpc_trans_plot(res, mpciterations, path)

clf; 
figure(1);
tcursor = 1:mpciterations;
subplot(3,1,1); plot(tcursor, res.x(tcursor,1), 'b-'); hold on;
plot(tcursor, path.xr(1:mpciterations), 'k--');
subplot(3,1,2); plot(tcursor, res.x(tcursor,3), 'b-'); hold on;
plot(tcursor, path.yr(1:mpciterations), 'k--');
subplot(3,1,3); plot(tcursor, res.x(tcursor,5), 'b-'); hold on;
plot(tcursor, path.zr(1:mpciterations), 'k--');
legend('Controlled', 'Reference');

figure(2);
subplot(2,1,1); plot(tcursor, res.u(tcursor,1),'k-'); 
legend('u_1');
subplot(2,1,2); plot(tcursor, res.u(tcursor,2),'b-'); hold on;
plot(tcursor, res.u(tcursor,3),'r-');
legend('u_x', 'u_y');

figure(3); 
plot3(path.xr(tcursor), ...
    path.yr(tcursor), path.zr(tcursor), 'k--'); hold on;
plot3(res.x(tcursor,1), res.x(tcursor,3), res.x(tcursor,5), 'b-');
legend('Controlled', 'Reference');
grid on;

figure(4);
subplot(3,1,1); plot(path.xr(tcursor)'-res.x(tcursor,1), 'k:'); hold on;
legend('\Delta_x');
subplot(3,1,2); plot(path.yr(tcursor)'-res.x(tcursor,3), 'b:'); hold on;
legend('\Delta_y');
subplot(3,1,3); plot(path.zr(tcursor)'-res.x(tcursor,5), 'r:');
legend('\Delta_z');

return