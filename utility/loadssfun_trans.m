function ss = loadssfun_trans()

phi = pi/4;
theta = pi/4;
m = 0.74;
U1 = 1;

ss.Axy = [0 1 0 0;...
    0 0 0 0;...
    0 0 0 1;...
    0 0 0 0];
ss.Bxy = [0 0; U1/m 0; 0 0; 0 U1/m];
ss.Qxy = diag([27000 1 21000 1]);
ss.Rxy = diag([1 1]);

ss.Az = [0 1;0 0];
ss.Bz = [0; 1/m*cos(phi)*cos(theta)];
ss.Qz = diag([21000 1]);
ss.Rz = diag(1);
end