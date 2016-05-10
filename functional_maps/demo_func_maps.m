function [Data_out] = demo_func_maps(Data, Para)
% Para.dimBasis = 1
if 1
    Data.shapes_ori = fm_orient_shapes(Data.shapes, 4);
end
if 1
    % Compute functional basis
    numShapes = length(Data.shapes);
    for shapeId = 1 : numShapes
        Data.basis{shapeId} = fm_man_made_basis(...
            Data.shapes{shapeId}, Para.dimBasis);
        fprintf('Finished computing basis for shape_%d\n', shapeId);
    end
end
if 1
    Data.shapeGraph = fm_shape_graph(Data.shapes_ori, Para);
end
if 1
    % Pairwise matching
    Data.input_pair_maps = fm_batch_pairwise_matching(...
        Data.shapes_ori, Data.shapeGraph, Para);
end
if 1
    % compute consistent functional maps
    Data.consistent_fmaps = fm_consistent_fmaps(Data, Para);
end
consistent_segmentations(Data, Para, Para.Camera);
Data_out = Data;