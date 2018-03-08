MATLAB toolbox of nonlinear model predictive control
   nmpc.m   nonlinear mpc main function from~\CITE{B-Grune-NMPC-2011}
   dataset  system prediction results using cgp <Based on the convloved gp toolbox>
   demos
   utils
   demo_sys_nmpc_quad    demo results by eps figures.

====================================================================================

NMPC based Nonlinear/Linear quadrotor system control

control strategy:  minmax/+terminal cost
    min_{u}{ max_{wn}{J(xn,un)}} = sum_n^n+N{(x-ref)*Q*(x-ref)' + (u-uref)*R*(u-uref)'}
    terminal cost = F(xn) = xN*P*xN';

 (demo_nmpc_quad_ss) -> linear 
    translational system consists of XY subsys and Z subsys;
   
 (demo_nmpc_quad)    -> nonlinear