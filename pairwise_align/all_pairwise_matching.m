function [SAMPLE, pairMatches] = all_pairwise_matching(Shapes, G_knn, Para_align)
%
% G_knn = knn_graph(Shapes, knn);

[rows, cols, vals] = find(G_knn);
ids = find(rows ~= cols);
rows = rows(ids)';
cols = cols(ids)';

for id = 1 : length(Shapes)
    sourcePoints_ori = sp_mesh_sampling(Shapes{id}, Para_align.numSamples);
    SAMPLE{id} = sourcePoints_ori;
    fprintf('Finished sampling %d.\n', id);
end

for pairId = 1 : length(rows)
    sId = rows(pairId);
    tId = cols(pairId);
    [ffd_opt, medianDis] = man_made_ffd_align(Shapes{sId}, Shapes{tId},...
        Para_align);
    
    % Compute close-point pairs
    C_source_ori = sp_ffd_basis_coeff(ffd_opt, SAMPLE{sId});
    sourcePoints = ffd_opt.ctrlPos_cur*C_source_ori';

    % Sample the target shape
    [IDX_s_t, dis_s_t] = knnsearch(SAMPLE{tId}', sourcePoints');
    sigma = median(dis_s_t);
    ids = find(dis_s_t < sigma*2);
    corres = [ids'; IDX_s_t(ids)'];
    pairMatches{pairId}.sId = sId;
    pairMatches{pairId}.tId = tId;
    pairMatches{pairId}.corres = corres;
    fprintf('Finished matching shape %d and shape %d.\n', sId, tId);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
