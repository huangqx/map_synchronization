function [G] = knn_graph(Shapes, knn)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute shape descriptors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Box.dimX = 10;
Box.dimY = 10;
Box.dimZ = 10;
Box.gridRes = 0.2;
Box.lowerCorner = [-1,-1,-1]';

dim = Box.dimX*Box.dimY*Box.dimZ;
numShapes = length(Shapes);
desscriptors = zeros(dim, numShapes);

for id = 1:numShapes
    Shape = Shapes{id};
    boxDim = max(Shape.vertexPoss')' - min(Shape.vertexPoss')';
    scaleX = 2/boxDim(1);
    scaleY = 2/boxDim(2);
    scaleZ = 2/boxDim(3);
    Shape.vertexPoss(1,:) = Shape.vertexPoss(1,:)*scaleX;
    Shape.vertexPoss(2,:) = Shape.vertexPoss(2,:)*scaleY;
    Shape.vertexPoss(3,:) = Shape.vertexPoss(3,:)*scaleZ;
    dess = lim_shape_descriptor(Shape, Box);
    desscriptors(:, id) = dess;
end

dim = numShapes;
Dis = zeros(dim, dim);
for i = 1:dim
    d = desscriptors(:,i)*ones(1,dim) - desscriptors;
    d = sqrt(sum(d.*d));
    Dis(i,:) = d;
end

G = sparse(dim, dim);

for i = 1:dim
    [s, ids] = sort(Dis(i,:));
    G(i, ids(2:(knn+1))) = 1;
end
G = max(G, G');

function [descriptors] = lim_shape_descriptor(Shape, Box)
%
points = lim_mesh_sampling(Shape, 16384);
cellCenters = zeros(3, Box.dimX*Box.dimY*Box.dimZ);
for i = 1:Box.dimX
    for j = 1:Box.dimY
        for k = 1:Box.dimZ
            id = (i-1)*Box.dimY*Box.dimZ + (j-1)*Box.dimZ + k;
            cellCenters(1, id) = Box.lowerCorner(1) + (i-0.5)*Box.gridRes;
            cellCenters(2, id) = Box.lowerCorner(2) + (j-0.5)*Box.gridRes;
            cellCenters(3, id) = Box.lowerCorner(3) + (k-0.5)*Box.gridRes;
        end
    end
end
[IDX, descriptors] = knnsearch(points', cellCenters');

function [points] = lim_mesh_sampling(Shape, numSamples)
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