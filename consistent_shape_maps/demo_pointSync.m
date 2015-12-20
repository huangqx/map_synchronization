%Load the data
[Data] = load_dataset('data\Fourleg\');

%Default parameters
Para.m0 = 24; % number of matched points
Para.m = 96; % number of samples where matched points are selected
Para.rootId = 1;

% Parameters used for solving the SDP
Para.mu_init = 0.1;
Para.nIterations = 800;
Para.rho = 1.015;

% Used in converting dense correspondences
Para.lambda = 1e-1;

% Compute samples
Data.SAMPLE = fps_sampling_all(Data, Para.m);

% Create SDP problem
[SDP] = csm_sdp_setup(Data, Para);

% Solve the SDP problem
[X] = csm_sdp_opt(SDP, Para);

% Obtained matched points
[consistentVertexIds] = csm_sdp_rounding(Data, X, Para);

% Interpolate into dense correspondences using functional maps
for i = 1:length(Data.shapes)
  Data.basis{i} = cotangent_basis(Data.shapes{i}, 32);
end
Data.opt_maps = batch_feature_2_vertex(Data, consistentVertexIds,...
    Para.lambda);

% Evaluate maps
[feaDisMat] = eval_precompt(Data);
curve_opt = eval_maps(Data.opt_maps, feaDisMat);
curve_initial = eval_maps(Data.initial_maps, feaDisMat);

plot(curve_initial(1:64), 'r');
hold on;
plot(curve_opt(1:64), 'b');

