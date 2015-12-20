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


% Perform sampling on the input shapes
SAMPLE = cell(1, length(Shapes));
for id = 1 : length(Shapes)
    SAMPLE{id} = sp_mesh_sampling(Shapes{id}, Para_align.numSamples);
    fprintf('Finished sampling %d.\n', id);
end

% Perform pair-wise affine matching between pairs of shapes
[unaryTerm, binaryTerm, edges, numOfStates, rootId] =...
    man_made_all_pairwise_matching(SAMPLE, G_knn, Para_align);

% Perform map synchronization to jointly optimize the transformations
% The current implementation uses the MRF formulation
% default parameters for running the MRF
paras = [400, 1e-9, 0.25];
sol = trws(numOfStates, unaryTerm, edges(1,:)-1, edges(2,:)-1,...
    binaryTerm, paras);
sol = [sol(1:(rootId-1)), 1, sol(rootId:length(sol))];
TP = sort(sol)

Shapes_aff = Shapes;
for id = 1 : length(Shapes)
    theta = (sol(id)-1)*2*pi/4/Para_align.numRotSamples;
    R = [cos(theta), 0, -sin(theta);
        0, 1, 0;
        sin(theta), 0, cos(theta)];
    Shapes_aff{id}.vertexPoss = R*Shapes_aff{id}.vertexPoss;
end
