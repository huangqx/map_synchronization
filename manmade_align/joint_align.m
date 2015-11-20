function [FFD] = joint_align(Shapes, SAMPLE, PAIRMATCH, Para_align)
%
numShapes = length(SAMPLE);
numPairs = length(PAIRMATCH);

for pairId = 1 : numPairs
    nc = size(PAIRMATCH{pairId}.corres, 2);
    PAIRMATCH{pairId}.weights = ones(1, nc);
end

numC_all = 0;
for pairId = 1 : numPairs
    numC = size(PAIRMATCH{pairId}.corres, 2);
    numC_all = numC_all + numC;
end

for outIter = 1:10
    % Perform recomputation
    for shapeId = 1 : numShapes
        weights{shapeId} = zeros(1, size(SAMPLE{shapeId}, 2));
    end

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

    SAMPLE_def = SAMPLE;
    for shapeId = 1 : numShapes
        FFD{shapeId} = sp_ffd_init_sym(Shapes{shapeId}, Para_align.gridRes);
        Term{shapeId}.b = sp_ffd_basis_coeff(FFD{shapeId}, SAMPLE{shapeId});
        w = weights{shapeId};
        W = sparse(1:length(w), 1:length(w), w);
    %    A = Ds*W*Ds' + Para.lambda_first*eye(dimX) + Para.lambda_smooth*FFD.H_smooth;
    %    b = Ds*W*Pt' + Para.lambda_first*FFD{shapeId}.ctrlPos_ori';
        dimX = size(Term{shapeId}.b, 2);
        Term{shapeId}.A = Term{shapeId}.b'*W*Term{shapeId}.b...
            + Para_align.lambda_first*eye(dimX)...
            + Para_align.lambda_smooth*FFD{shapeId}.H_smooth;
        Term{shapeId}.b = Term{shapeId}.b';
    end

    for iter = 1:32
        % Perform free-form deformation to obtain deformed positions
        for shapeId = 1 : numShapes
            SAMPLE_def{shapeId} = FFD{shapeId}.ctrlPos_cur*Term{shapeId}.b;
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
                Para_align.lambda_first*FFD{shapeId}.ctrlPos_ori';
            FFD{shapeId}.ctrlPos_cur = (Term{shapeId}.A\b)';
        end
    end
    % Reweight the correspondences
    for shapeId = 1 : numShapes
        SAMPLE_def{shapeId} = FFD{shapeId}.ctrlPos_cur*Term{shapeId}.b;
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


% Compute the coefficient of 3D positions with respect to a free-form
% deformation structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C] = sp_ffd_basis_coeff(FFD, poss)
%
poss = double(poss);

tX = (poss(1,:) - FFD.lc(1))/FFD.gridRes;
idX = floor(tX) + 1;
tX = tX - (idX-1);

tY = (poss(2,:) - FFD.lc(2))/FFD.gridRes;
idY = floor(tY) + 1;
tY = tY - (idY-1);

tZ = (poss(3,:) - FFD.lc(3))/FFD.gridRes;
idZ = floor(tZ) + 1;
tZ = tZ - (idZ-1);

nP = size(poss, 2);

rows = ones(8,1)*(1:nP);
cols = zeros(8, nP);
vals = zeros(8, nP);

cols(1,:) = (idY-1)*FFD.dimX*FFD.dimZ + (idX-1)*FFD.dimZ + idZ;
vals(1,:) = (1-tX).*(1-tY).*(1-tZ);

cols(2,:) = (idY-1)*FFD.dimX*FFD.dimZ + (idX-1)*FFD.dimZ + idZ + 1;
vals(2,:) = (1-tX).*(1-tY).*tZ;

cols(3,:) = (idY)*FFD.dimX*FFD.dimZ + (idX-1)*FFD.dimZ + idZ;
vals(3,:) = (1-tX).*tY.*(1-tZ);

cols(4,:) = (idY)*FFD.dimX*FFD.dimZ + (idX-1)*FFD.dimZ + idZ + 1;
vals(4,:) = (1-tX).*tY.*tZ;

cols(5,:) = (idY-1)*FFD.dimX*FFD.dimZ + (idX)*FFD.dimZ + idZ;
vals(5,:) = tX.*(1-tY).*(1-tZ);

cols(6,:) = (idY-1)*FFD.dimX*FFD.dimZ + (idX)*FFD.dimZ + idZ + 1;
vals(6,:) = tX.*(1-tY).*tZ;

cols(7,:) = (idY)*FFD.dimX*FFD.dimZ + (idX)*FFD.dimZ + idZ;
vals(7,:) =  tX.*tY.*(1-tZ);

cols(8,:) = (idY)*FFD.dimX*FFD.dimZ + (idX)*FFD.dimZ + idZ + 1;
vals(8,:) =  tX.*tY.*tZ;

C = sparse(rows, cols, vals, nP, FFD.dimX*FFD.dimY*FFD.dimZ);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initilaize the free-from deformation structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FFD] = sp_ffd_init_sym(Shape, gridRes)
% Initialize the FFD data structure, given a grid resolution
% Input arguments:
%   Shape: the input shape
%   gridRes: the resolution of the grid
% Ouput argument:
%   FFD: the initialize grid-structure

% Compute the bounding box
boxMin = double(min(Shape.vertexPoss'));
boxMax = double(max(Shape.vertexPoss'));
boxCenter = (boxMin + boxMax)/2;
boxSize = boxMax - boxMin;

% Initialize the grid
dimX = double(floor(boxSize(1)/gridRes) + 2);
dimY = double(floor(boxSize(2)/gridRes) + 2);
dimZ = double(floor(boxSize(3)/gridRes) + 2);

FFD.lc(1) = boxCenter(1) - gridRes*(dimX-1)/2;
FFD.lc(2) = boxCenter(2) - gridRes*(dimY-1)/2;
FFD.lc(3) = boxCenter(3) - gridRes*(dimZ-1)/2;
FFD.gridRes = gridRes;
FFD.dimX = dimX;
FFD.dimY = dimY;
FFD.dimZ = dimZ;

% Initilize the gridRes solution
FFD.ctrlPos_ori = zeros(3, dimX*dimY*dimZ);
IDX = 1:(dimX*dimY*dimZ);

offsetsY = floor((IDX-1)/(dimX*dimZ)) + 1;
offsetsZ = IDX - (offsetsY-1)*dimX*dimZ;
offsetsX = floor((offsetsZ-1)/dimZ) + 1;
offsetsZ = offsetsZ - (offsetsX-1)*dimZ;

FFD.ctrlPos_ori(1,:) = FFD.lc(1) + (offsetsX-1)*gridRes;
FFD.ctrlPos_ori(2,:) = FFD.lc(2) + (offsetsY-1)*gridRes;
FFD.ctrlPos_ori(3,:) = FFD.lc(3) + (offsetsZ-1)*gridRes;
FFD.ctrlPos_cur = FFD.ctrlPos_ori;

FFD.offsetsX = offsetsX;
FFD.offsetsY = offsetsY;
FFD.offsetsZ = offsetsZ;

if mod(dimX, 2) == 0
    colids_1 = find(offsetsX < (dimX+1)/2);
    colids_2 = (offsetsY(colids_1)-1)*dimX*dimZ +...
        (dimX - offsetsX(colids_1))*dimZ +...
        offsetsZ(colids_1);
    
    rows = kron(1:length(colids_1), ones(1,2));
    cols = reshape([colids_1;colids_2], [1, 2*length(colids_1)]);
    vals = kron(ones(1,length(colids_1)), [1,1]);
    FFD.Ax = sparse(rows, cols, vals);
    FFD.bx = zeros(size(FFD.Ax,1),1);
    
    vals = kron(ones(1,length(colids_1)), [1,-1]);
    FFD.Az = sparse(rows, cols, vals);
    FFD.bz = zeros(size(FFD.Az,1),1);
    
    colids_1 = find(offsetsY == 1);
    colids_2 = find(offsetsY == dimY);
    rows = 1:(2*dimX*dimZ);
    cols = double([colids_1, colids_2]);
    vals = ones(1, 2*dimX*dimZ);
    Ay1 = sparse(rows, cols, vals, 2*dimX*dimZ, dimX*dimY*dimZ);
    yMin = FFD.lc(2);
    yMax = yMin + (dimY-1)*gridRes;
    by1 = [ones(1, length(colids_1))*yMin, ones(1, length(colids_2))*yMax]';
    
    colids_1 = find(offsetsX < (dimX+1)/2 & offsetsY > 1 & offsetsY < dimY);
    colids_2 = (offsetsY(colids_1)-1)*dimX*dimZ +...
    (dimX - offsetsX(colids_1))*dimZ +...
    offsetsZ(colids_1);
    rows = kron(1:length(colids_1), ones(1,2));
    cols = reshape([colids_1;colids_2], [1, 2*length(colids_1)]);
    vals = kron(ones(1,length(colids_1)), [1,-1]);
    Ay2 = sparse(rows, cols, vals, length(colids_1), dimX*dimY*dimZ);
    by2 = zeros(size(Ay2,1),1);
    
    FFD.Ay = [Ay1;Ay2];
    FFD.by = [by1;by2];
else
    cols = find(offsetsX == (dimX+1)/2);
    rows = 1:length(cols);
    vals = ones(1, length(cols));
    
    colids_1 = find(offsetsX < dimX/2);
    colids_2 = (offsetsY(colids_1)-1)*dimX*dimZ +...
        (dimX - offsetsX(colids_1))*dimZ +...
        offsetsZ(colids_1);
    
    rows = [rows, length(rows)+kron(1:length(colids_1), ones(1,2))];
    cols = [cols, reshape([colids_1;colids_2], [1, 2*length(colids_1)])];
    vals = [vals, kron(ones(1,length(colids_1)), [1,1])];
    FFD.Ax = sparse(rows, cols, vals);
    FFD.bx = zeros(size(FFD.Ax,1),1);
    
    rows = [kron(1:length(colids_1), ones(1,2))];
    cols = [reshape([colids_1;colids_2], [1, 2*length(colids_1)])];
    vals = kron(ones(1,length(colids_1)), [1,-1]);
    FFD.Az = sparse(rows, cols, vals);
    FFD.bz = zeros(size(FFD.Az,1),1);
    
    colids_1 = find(offsetsY == 1);
    colids_2 = find(offsetsY == dimY);
    rows = 1:(2*dimX*dimZ);
    cols = double([colids_1, colids_2]);
    vals = ones(1, 2*dimX*dimZ);
    Ay1 = sparse(rows, cols, vals, 2*dimX*dimZ, dimX*dimY*dimZ);
    yMin = FFD.lc(2);
    yMax = yMin + (dimY-1)*gridRes;
    by1 = [ones(1, length(colids_1))*yMin, ones(1, length(colids_2))*yMax]';
    
    colids_1 = find(offsetsX < dimX/2 & offsetsY > 1 & offsetsY < dimY);
    colids_2 = (offsetsY(colids_1)-1)*dimX*dimZ +...
    (dimX - offsetsX(colids_1))*dimZ +...
    offsetsZ(colids_1);
    rows = kron(1:length(colids_1), ones(1,2));
    cols = reshape([colids_1;colids_2], [1, 2*length(colids_1)]);
    vals = kron(ones(1,length(colids_1)), [1,-1]);
    Ay2 = sparse(rows, cols, vals, length(colids_1), dimX*dimY*dimZ);
    by2 = zeros(size(Ay2,1),1);
    
    FFD.Ay = [Ay1;Ay2];
    FFD.by = [by1;by2];
end

% The smoothness term
subset = find(offsetsX < dimX -1);
vIds0 = (offsetsY(subset)-1)*dimX*dimZ + (offsetsX(subset)-1)*dimZ + offsetsZ(subset);
vIds1 = (offsetsY(subset)-1)*dimX*dimZ + (offsetsX(subset))*dimZ + offsetsZ(subset);
vIds2 = (offsetsY(subset)-1)*dimX*dimZ + (offsetsX(subset)+1)*dimZ + offsetsZ(subset);

rows = kron(1:length(vIds0), ones(1,3));
cols = [vIds0;vIds1;vIds2];
vals = kron(ones(1,length(vIds0)), [1,-2,1]);
Jx = sparse(rows, cols, vals, length(vIds0), dimX*dimY*dimZ);

subset = find(offsetsY < dimY -1);
vIds0 = (offsetsY(subset)-1)*dimX*dimZ + (offsetsX(subset)-1)*dimZ + offsetsZ(subset);
vIds1 = (offsetsY(subset))*dimX*dimZ + (offsetsX(subset)-1)*dimZ + offsetsZ(subset);
vIds2 = (offsetsY(subset)+1)*dimX*dimZ + (offsetsX(subset)-1)*dimZ + offsetsZ(subset);

rows = kron(1:length(vIds0), ones(1,3));
cols = [vIds0;vIds1;vIds2];
vals = kron(ones(1,length(vIds0)), [1,-2,1]);
Jy = sparse(rows, cols, vals, length(vIds0), dimX*dimY*dimZ);

subset = find(offsetsZ < dimZ -1);
vIds0 = (offsetsY(subset)-1)*dimX*dimZ + (offsetsX(subset)-1)*dimZ + offsetsZ(subset);
vIds1 = (offsetsY(subset)-1)*dimX*dimZ + (offsetsX(subset)-1)*dimZ + offsetsZ(subset)+1;
vIds2 = (offsetsY(subset)-1)*dimX*dimZ + (offsetsX(subset)-1)*dimZ + offsetsZ(subset)+2;

rows = kron(1:length(vIds0), ones(1,3));
cols = [vIds0;vIds1;vIds2];
vals = kron(ones(1,length(vIds0)), [1,-2,1]);
Jz = sparse(rows, cols, vals, length(vIds0), dimX*dimY*dimZ);

J_smooth = [Jx;Jy;Jz];
FFD.H_smooth = J_smooth'*J_smooth;
