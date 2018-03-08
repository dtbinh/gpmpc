% time varying u=Kx mpc

clc; clear; close all;
pathtype = 'Elliptical';    % Lorenz   Elliptical   Helix   crown   spraypool
mpciterations = 199;   % 100 120  150  -> H=10/20

[gppred, refpred, varargout] = ...
    pgp_tvKmpc_quad(pathtype, mpciterations, 'TRAN');
