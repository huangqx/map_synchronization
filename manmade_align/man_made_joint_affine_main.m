function [Shapes_aff] = man_made_joint_affine_main(Shapes, Para_align)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimize a free-from deformation for each shape so that the input shapes
% are aligned in a world coordinate system
% Input arguments:
%       Shapes:     the input shapes
%       Para_align: the parameters used in pair-wise alignment
% Output argument:
%       Shapes_aff: the transformed shapes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Computing the knn shape graph...\n');
% Compute a shape-graph that aligns all the shapes
G_knn = man_made_knn_graph(Shapes, Para_align.knn);

% Perform pair-wise affine matching between pairs of shapes
[SAMPLE, PAIRMATCH] = man_made_all_pairwise_match(Shapes, G_knn, Para_align);

