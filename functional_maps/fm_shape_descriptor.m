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