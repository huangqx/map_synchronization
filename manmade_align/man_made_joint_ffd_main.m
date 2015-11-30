function [FFDs_opt] = man_made_joint_ffd_main(Shapes, Para_align)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimize a free-from deformation for each shape so that the input shapes
% are aligned in a world coordinate system
% Input arguments:
%       Shapes:     the input shapes
%       Para_align: the parameters used in pair-wise alignment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Computing the knn shape graph...\n');
G_knn = man_made_knn_graph(Shapes, Para_align.knn);

% Perform pair-wise matching between all pairs of shapes
fprintf('Performing all pair-wise matching...\n');
[SAMPLE, PAIRMATCH] = man_made_all_pairwise_align(Shapes, G_knn, Para_align);

% Perform joint-alignment
fprintf('Performing joint alignment...\n');
FFDs_opt = man_made_joint_align(Shapes, SAMPLE, PAIRMATCH, Para_align);
