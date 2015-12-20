function [Para_csm] = set_csm_parameters()

Para_csm.m0 = 24; % number of matched points
Para_csm.m = 96; % number of samples where matched points are selected
Para_csm.rootId = 1;

% Parameters used for solving the SDP
Para_csm.mu_init = 0.1;
Para_csm.nIterations = 800;
Para_csm.rho = 1.015;

% Used in converting dense correspondences
Para_csm.lambda = 1e-1;
