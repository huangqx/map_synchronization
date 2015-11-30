function [SAMPLE, PAIRMATCH] = man_made_all_pairwise_align(Shapes,...% The input shapes
    G_knn,... % The shape graph
    Para_align) % The parameters
% This function performs all pair-wise matching along pairs of shapes
% specified by a shape graph

% Extract pairs of shapes to be aligned
[rows, cols, vals] = find(G_knn);
ids = find(rows ~= cols);
rows = rows(ids)';
cols = cols(ids)';

% Perform sampling on the input shapes
for id = 1 : length(Shapes)
    SAMPLE{id} = sp_mesh_sampling(Shapes{id}, Para_align.numSamples);
    FFDs{id} = sp_ffd_init_sym(Shapes{id}, Para_align.gridRes);
    fprintf('Finished sampling %d.\n', id);
end

% Perform pair-wise alignment
for pairId = 1 : length(rows)
    PAIRMATCH{pairId} = pairwise_matching(...
        FFDs{pairId},...
        SAMPLE{rows(pairId)},...
        SAMPLE{cols(pairId)},...
        rows(pairId),...
        cols(pairId),...
        Para_align);
    fprintf('Finished matching shape %d and shape %d.\n',  rows(pairId),  cols(pairId));
end

function [match] = pairwise_matching(...
    sourceFFD,...
    sourceSample,...
    targetSample,...
    sId, tId,...
    Para_align)
%
[ffd_opt, medianDis] = man_made_pairwise_ffd_align(...
    sourceFFD,...
    sourceSample,...
    targetSample,...
    Para_align);
    
% After alignment, compute the closest point pairs
C_source_ori = sp_ffd_basis_coeff(ffd_opt, sourceSample); % FFD coefficient
sourcePoints = ffd_opt.ctrlPos_cur*C_source_ori'; % Current point positions
[IDX_s_t, dis_s_t] = knnsearch(targetSample', sourcePoints'); % Compute nearest neighbors
sigma = median(dis_s_t); % Use robust thresholding to find correspondences
ids = find(dis_s_t < sigma*2);
corres = [ids'; IDX_s_t(ids)'];
    
% Store the correspondences
match.sId = sId;
match.tId = tId;
match.corres = corres;




