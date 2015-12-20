function [opt_maps] = csm_main_func(Data, Para_csm)
% This is the main function for optimizing consistent shape maps. The
% function takes a collection of shapes and maps computed between pairs of
% shapes as input and output optimized maps.

% Compute samples
Data.SAMPLE = fps_sampling_all(Data, Para_csm.m);

% Create SDP problem
[SDP] = csm_sdp_setup(Data, Para_csm);

% Solve the SDP problem
[X] = csm_sdp_opt(SDP, Para_csm);

% Obtained matched points
[consistentVertexIds] = csm_sdp_rounding(Data, X, Para_csm);

% Interpolate into dense correspondences using functional maps
for i = 1:length(Data.shapes)
  Data.basis{i} = cotangent_basis(Data.shapes{i}, 32);
end

opt_maps = batch_feature_2_vertex(Data, consistentVertexIds,...
    Para_csm.lambda);

% Evaluate maps
[feaDisMat] = eval_precompt(Data);
curve_opt = eval_maps(opt_maps, feaDisMat);
curve_initial = eval_maps(Data.initial_maps, feaDisMat);

plot(curve_initial(1:64), 'r');
hold on;
plot(curve_opt(1:64), 'b');

