function [FFDs] = man_made_joint_align(Shapes, SAMPLE, PAIRMATCH, Para_align)
% Given the correspondences computed between all pairs of shapes
% Align all the shapes in the world coordinate system
% Input arguments:
%       Shapes: the input shapes
%       SAMPLE: the samples placed on each shape
%       PAIRMATCH: the pre-computed correspondences between pairs of shapes
%       Para_align: the parameters
% Output argument:
%        FFDs{shapeId}: the optimized free-form deformation of each shape

numShapes = length(SAMPLE);
numPairs = length(PAIRMATCH);

% Allocate space to store the pair-wise correspondence weights
for pairId = 1 : numPairs
    nc = size(PAIRMATCH{pairId}.corres, 2);
    PAIRMATCH{pairId}.weights = ones(1, nc);
end

% Count the number of total correspondences
numC_all = 0;
for pairId = 1 : numPairs
    numC = size(PAIRMATCH{pairId}.corres, 2);
    numC_all = numC_all + numC;
end

for outIter = 1:Para_align.numIterations_outer
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Precompute the number of correspondences acted on each sample
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Allocate space for cumulative weights
    for shapeId = 1 : numShapes
        weights{shapeId} = zeros(1, size(SAMPLE{shapeId}, 2));
    end

    % Pre-compute the weights
    for pairId = 1 : numPairs
        nc = size(PAIRMATCH{pairId}.corres, 2);
        PAIRMATCH{pairId}.midPoints = zeros(3, nc);
        sId = PAIRMATCH{pairId}.sId;
        tId = PAIRMATCH{pairId}.tId;
        ids1 = PAIRMATCH{pairId}.corres(1,:);
        ids2 = PAIRMATCH{pairId}.corres(2,:);
        weights{sId}(ids1) = weights{sId}(ids1) + PAIRMATCH{pairId}.weights;
        weights{tId}(ids2) = weights{tId}(ids2) + PAIRMATCH{pairId}.weights;
    end

    % Pre-compute the deformation structure used for each shape
    SAMPLE_def = SAMPLE;
    for shapeId = 1 : numShapes
        FFDs{shapeId} = sp_ffd_init_sym(Shapes{shapeId}, Para_align.gridRes);
        Term{shapeId}.b = sp_ffd_basis_coeff(FFDs{shapeId}, SAMPLE{shapeId});
        w = weights{shapeId};
        W = sparse(1:length(w), 1:length(w), w);
        dimX = size(Term{shapeId}.b, 2);
        Term{shapeId}.A = Term{shapeId}.b'*W*Term{shapeId}.b...
            + Para_align.lambda_first*eye(dimX)...
            + Para_align.lambda_smooth*FFDs{shapeId}.H_smooth;
        Term{shapeId}.b = Term{shapeId}.b';
    end
    
    % Peform alternating optimization to optimize the deformation on each
    % shape
    for iter = 1:Para_align.numIterations_alternate
        % Compute the deformed positions of sample points
        for shapeId = 1 : numShapes
            SAMPLE_def{shapeId} = FFDs{shapeId}.ctrlPos_cur*Term{shapeId}.b;
        end
    
        % Compute the intermediate points
        sqrDis = 0;
        for pairId = 1 : numPairs
            sId = PAIRMATCH{pairId}.sId;
            tId = PAIRMATCH{pairId}.tId;
            poss_s = SAMPLE_def{sId}(:, PAIRMATCH{pairId}.corres(1,:));
            poss_t = SAMPLE_def{tId}(:, PAIRMATCH{pairId}.corres(2,:));
            PAIRMATCH{pairId}.midPoints = (poss_s + poss_t)/2;
            d = poss_s - poss_t;
            sqrDis = sqrDis + sum(sum(d.*d).*PAIRMATCH{pairId}.weights);
        end
        fprintf('sqrDis = %f.\n', sqrDis);
    
        % Perform the alignment
        for shapeId = 1 : numShapes
            TP{shapeId} = zeros(3, size(SAMPLE{shapeId}, 2));
        end
        for pairId = 1 : numPairs
            sId = PAIRMATCH{pairId}.sId;
            tId = PAIRMATCH{pairId}.tId;
            ids1 = PAIRMATCH{pairId}.corres(1,:);
            ids2 = PAIRMATCH{pairId}.corres(2,:);
            buf = PAIRMATCH{pairId}.midPoints.*(ones(3,1)*PAIRMATCH{pairId}.weights);
            TP{sId}(:, ids1) = TP{sId}(:, ids1) + buf;
            TP{tId}(:, ids2) = TP{tId}(:, ids2) + buf;
        end
        for shapeId = 1 : numShapes
            b = Term{shapeId}.b*TP{shapeId}' +...
                Para_align.lambda_first*FFDs{shapeId}.ctrlPos_ori';
            FFDs{shapeId}.ctrlPos_cur = (Term{shapeId}.A\b)';
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reweight the correspondences
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for shapeId = 1 : numShapes
        SAMPLE_def{shapeId} = FFDs{shapeId}.ctrlPos_cur*Term{shapeId}.b;
    end
    %
    sqrDisVec = zeros(1, numC_all);
    off = 0;
    for pairId = 1 : numPairs
        sId = PAIRMATCH{pairId}.sId;
        tId = PAIRMATCH{pairId}.tId;
        poss_s = SAMPLE_def{sId}(:, PAIRMATCH{pairId}.corres(1,:));
        poss_t = SAMPLE_def{tId}(:, PAIRMATCH{pairId}.corres(2,:));
        d = poss_s - poss_t;
        d = sum(d.*d);
        nc = length(d);
        sqrDisVec((off+1):(off+nc)) = d;
        off = off + nc;
    end
    sigma = median(sqrDisVec);
    for pairId = 1 : numPairs
        sId = PAIRMATCH{pairId}.sId;
        tId = PAIRMATCH{pairId}.tId;
        poss_s = SAMPLE_def{sId}(:, PAIRMATCH{pairId}.corres(1,:));
        poss_t = SAMPLE_def{tId}(:, PAIRMATCH{pairId}.corres(2,:));
        d = poss_s - poss_t;
        d = sum(d.*d);
        PAIRMATCH{pairId}.weights = sqrt(sigma)./sqrt(sigma + d);
    end
end

