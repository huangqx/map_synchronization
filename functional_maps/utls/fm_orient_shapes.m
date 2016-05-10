function [Shapes_cons] = fm_orient_shapes(Shapes, numCandOrients)
%
gridDim = 10;
%
numShapes = length(Shapes);
dim = numShapes*numCandOrients;
dess = zeros(gridDim*gridDim*gridDim, dim);
for shapeId = 1 : numShapes
    for id = 1 : numCandOrients
        theta = 2*pi*(id-1)/numCandOrients;
        R = [cos(theta), 0, -sin(theta); 0, 1, 0; sin(theta), 0, cos(theta)];
        Shape = Shapes{shapeId};
        Shape.vertexPoss = R*Shape.vertexPoss;
        dess(:, (shapeId-1)*numCandOrients + id) =...
            fm_shape_descriptor(Shape, gridDim);
    end
end
%
dissScores = zeros(dim, dim);
for id = 1 : dim
    for id2 = (id+1) : dim
        dissScores(id, id2) = norm(dess(:, id) - dess(:, id2));
        dissScores(id2, id) = dissScores(id, id2);
    end
end
%
diss = zeros(numShapes, numShapes);
for id = 1:numShapes
    ids = ((id-1)*numCandOrients + 1):(id*numCandOrients);
    for id2 = 1:numShapes
        ids2 = ((id2-1)*numCandOrients + 1):(id2*numCandOrients);
        diss(id, id2) = min(min(dissScores(ids, ids2)));
    end
end
diss = reshape(diss, [1, numShapes*numShapes]);
sigma = median(diss);
scores = exp(-dissScores.*dissScores/2/sigma/sigma);
IDX = (numCandOrients+1):dim;
b = scores(IDX, 1)*numShapes;
A = scores(IDX, IDX);
for i = 1:(numShapes-1)
    ids2 = ((i-1)*numCandOrients+1):(i*numCandOrients);
    A(ids2, ids2) = 0;
end
%
x = ones(length(IDX), 1)/sqrt(numCandOrients);
for iter = 1 : 32
    x = A*x + b;
    x = reshape(x, [numCandOrients, (numShapes-1)]);
    n = sqrt(sum(x.*x));
    x = x./(ones(numCandOrients,1)*n);
    x = reshape(x, [length(IDX), 1]);
end
x = reshape(x, [numCandOrients, numShapes-1]);
[s, orientIds] = max(x);

for shapeId = 2 : numShapes
    theta = 2*pi*(orientIds(shapeId-1)-1)/numCandOrients;
    R = [cos(theta), 0, -sin(theta); 0, 1, 0; sin(theta), 0, cos(theta)];
    Shape = Shapes{shapeId};
    Shape.vertexPoss = R*Shape.vertexPoss;
    Shapes{shapeId} = Shape;
end
Shapes_cons = Shapes;