function [SAMPLE, PAIRMATCH] = man_made_all_pairwise_matching(...
    Shapes,...% The input shapes
    G_knn,... % The shape graph
    Para_align) % The parameters
%
% This function performs all pair-wise affine matching along pairs of shapes
% specified by a shape graph

% Extract pairs of shapes to be aligned
[rows, cols, vals] = find(G_knn);
ids = find(rows ~= cols);
rows = rows(ids)';
cols = cols(ids)';

% Perform sampling on the input shapes
for id = 1 : length(Shapes)
    SAMPLE{id} = sp_mesh_sampling(Shapes{id}, Para_align.numSamples);
    fprintf('Finished sampling %d.\n', id);
end

% Perform pair-wise matching
for pairId = 1 : length(rows)
    sId = rows(pairId);
    tId = cols(pairId);
    aff = man_made_pairwise_affine_matching(...
        SAMPLE{sId},...
        SAMPLE{tId},...
        Para_align);
    
    PAIRMATCH{pairId}.sId = sId;
    PAIRMATCH{pairId}.tId = tId;
    PAIRMATCH{pairId}.aff_trans = aff;
    fprintf('Finished matching shape %d and shape %d.\n', sId, tId);
end

% Perform map synchronization to jointly optimize the transformations
% The current implementation uses the MRF formulation 

