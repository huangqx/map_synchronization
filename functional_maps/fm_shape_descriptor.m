function [dess] = fm_shape_descriptor(Shape, gridDim)
% compute the volumetric shape descriptor with grid dimension 'gridDim'
bbox = max(Shape.vertexPoss')' - min(Shape.vertexPoss')';
center = (max(Shape.vertexPoss')' + min(Shape.vertexPoss')')/2;
%
numV = size(Shape.vertexPoss, 2);
Shape.vertexPoss = Shape.vertexPoss - center*ones(1, numV);
for i = 1:3
    Shape.vertexPoss(i, :) = Shape.vertexPoss(i, :)/bbox(i);
end
%
cellCenters = zeros(3, gridDim*gridDim*gridDim);
for i = 1:dim
    for j = 1:dim
        for k = 1:dim
            cellId = (i-1)*dim*dim + (j-1)*dim + k;
            cellCenters(1, cellId) = (i-0.5)/dim;
            cellCenters(2, cellId) = (j-0.5)/dim;
            cellCenters(3, cellId) = (k-0.5)/dim;
        end
    end
end
% Perform uniform sampling of the shape
numSamples = 32768;
points = fm_mesh_sampling(Shape, numSamples);
%
dim = gridDim*gridDim*gridDim;
dess = zeros(dim, 1);
%
d = kron(cellCenters, ones(1, numSamples)) - kron(ones(1, dim), points);
d = sum(d.*d);
d = reshape(d, [numSamples, dim]);
dess = sqrt(min(d));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform uniform sampling of each mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [points] = fm_mesh_sampling(Shape, numSamples)
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
