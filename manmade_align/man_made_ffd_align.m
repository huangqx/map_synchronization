function [ffd_opt, medianDis] = man_made_ffd_align(...
    sourceShape, targetShape, Para)
% Optimize a free-form deformation that aligns the source shape to the
% target shape

% Initialize the free-form deformation structure
FFD = sp_ffd_init_sym(sourceShape, Para.gridRes);

% Sample the source shape
sourcePoints_ori = sp_mesh_sampling(sourceShape, Para.numSamples);

% Sample the target shape
targetPoints = sp_mesh_sampling(targetShape, Para.numSamples);

% The FFD coefficient of the source points
C_source_ori = sp_ffd_basis_coeff(FFD, sourcePoints_ori);

% Applies non-rigid 
for iter = 1:Para.numIterations
    sourcePoints = FFD.ctrlPos_cur*C_source_ori';
    % Compute correspondences
    [Corres, medianDis] = sp_closest_point_corres(sourcePoints,...
        targetPoints);
    
    % Deform the shape accordingly
    nC = size(Corres, 2);
    W = sparse(1:nC, 1:nC, Corres(3,:));
    
    Ds = C_source_ori(Corres(1,:),:)';
    Pt = targetPoints(:, Corres(2,:));
    
    dimX = size(Ds, 1);
    A = Ds*W*Ds' + Para.lambda_first*eye(dimX) + Para.lambda_smooth*FFD.H_smooth;
    b = Ds*W*Pt' + Para.lambda_first*FFD.ctrlPos_ori';
    FFD.ctrlPos_cur = (A\b)';
end
%
ffd_opt = FFD;

sourceShape_def = sourceShape;
C_source_ori = sp_ffd_basis_coeff(FFD, sourceShape.vertexPoss);
sourceShape_def.vertexPoss = FFD.ctrlPos_cur*C_source_ori';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute closest point pairs between three point clouds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Corres, medianDis] = sp_closest_point_corres(sourcePoints,...
    targetPoints)
%
% Compute cloest point pairs in both directions
[footIds_s_in_t, d_s_in_t] = knnsearch(targetPoints', sourcePoints');
[footIds_t_in_s, d_t_in_s] = knnsearch(sourcePoints', targetPoints');

dis = [d_s_in_t', d_t_in_s'];
Corres_s = [1:length(footIds_s_in_t), footIds_t_in_s'];
Corres_t = [footIds_s_in_t', 1:length(footIds_t_in_s)];

Corres = [Corres_s; Corres_t];
d = sort(dis);
sigma = d(floor(length(d)*0.85));
sigma = max(5e-2, sigma);
weights = sigma*sigma./(sigma*sigma + dis.*dis);

Corres = [Corres; weights];
medianDis = sqrt(mean(dis.*dis));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform uniform sampling of each mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [points] = sp_mesh_sampling(Shape, numSamples)
%
numFaces = size(Shape.faceVIds, 2);
faceAreas = zeros(1, numFaces);
for faceId = 1:numFaces
    p1 = Shape.vertexPoss(:, Shape.faceVIds(1, faceId));
    p2 = Shape.vertexPoss(:, Shape.faceVIds(2, faceId));
    p3 = Shape.vertexPoss(:, Shape.faceVIds(3, faceId));
    e12 = p1 - p2;
    e13 = p1 - p3;
    e23 = p2 - p3;
    a = e12'*e12;
    b = e13'*e13;
    c = e23'*e23;
    faceAreas(faceId) = sqrt(2*(a*b+a*c+b*c) - (a*a+b*b+c*c))/4;
    if faceId > 1
        faceAreas(faceId) = faceAreas(faceId) + faceAreas(faceId-1);
    end
end

faceAreas = faceAreas/faceAreas(numFaces);
sample_ts = sort(rand(1, numSamples));
points = zeros(3, numSamples);

faceId = 1;
for sId = 1 : numSamples
    while sample_ts(sId) > faceAreas(faceId)
        faceId = faceId + 1;
    end
    p1 = Shape.vertexPoss(:, Shape.faceVIds(1,faceId));
    p2 = Shape.vertexPoss(:, Shape.faceVIds(2,faceId));
    p3 = Shape.vertexPoss(:, Shape.faceVIds(3,faceId));
    r1 = rand(1,1);
    r2 = rand(1,1);
    t1 = (1-sqrt(r1));
    t2 = sqrt(r1)*(1-r2);
    t3 = sqrt(r1)*r2;
    points(:, sId) = t1*p1 + t2*p2 + t3*p3;
end