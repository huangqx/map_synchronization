function [input_pair_maps] = fm_batch_pairwise_matching(Shapes, shapeGraph, Para)
% samples are given by centers of face triangles
numShapes = length(Shapes);
for shapeId = 1 : numShapes
    Shape = Shapes{shapeId};
    points1 = Shape.vertexPoss(:, Shape.faceVIds(1,:));
    points2 = Shape.vertexPoss(:, Shape.faceVIds(2,:));
    points3 = Shape.vertexPoss(:, Shape.faceVIds(3,:));
    SAMPLE{shapeId} = (points1 + points2 + points3)/3;
end
% intialize the FFDS
for id = 1 : length(Shapes)
    FFDs{id} = fm_ffd_init_sym(Shapes{id}, Para.gridRes);
end
% Extract edges
[rows, cols, vals] = find(shapeGraph);
% Perform pair-wise alignment
input_pair_maps = cell(1, length(rows));
parfor pairId = 1 : length(rows)
    input_pair_maps{pairId} = pairwise_matching(...
        FFDs{rows(pairId)},...
        SAMPLE{rows(pairId)},...
        SAMPLE{cols(pairId)},...
        rows(pairId),...
        cols(pairId),...
        Para);
    fprintf('Finished matching shape %d and shape %d.\n',  rows(pairId),  cols(pairId));
end

function [match] = pairwise_matching(...
    sourceFFD,...
    sourceSample,...
    targetSample,...
    sId, tId,...
    Para)
%
[ffd_opt, medianDis] = fm_pairwise_ffd_align(...
    sourceFFD,...
    sourceSample,...
    targetSample,...
    Para);
    
% After alignment, compute the closest point pairs
C_source_ori = fm_ffd_basis_coeff(ffd_opt, sourceSample); % FFD coefficient
sourcePoints = ffd_opt.ctrlPos_cur*C_source_ori'; % Current point positions
[IDX_s_t, dis_s_t] = knnsearch(targetSample', sourcePoints'); % Compute nearest neighbors
sigma = median(dis_s_t); % Use robust thresholding to find correspondences
ids = find(dis_s_t < sigma*3);
dis_s_t = dis_s_t(ids)';
corres = [ids'; IDX_s_t(ids)'; exp(-dis_s_t.^2/2/sigma/sigma)];
    
% Store the correspondences
match.sId = sId;
match.tId = tId;
match.corres = corres;