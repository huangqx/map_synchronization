function [Basis] = manmade_basis(Shape, dimBasis, sampleDensity)
% Compute the functional basis of a triangular mesh Shape
% Shape.vertexPoss and Shape.faceVIds
% Note that the Shape may contain many disconnected components and the
% triangular faces maybe irregular
% Other arguments:
%   dimBasis: the dimension of the functional space
%   sampleDensity: the density of the additional samples 
%
Samples = mesh_sampling(Shape, sampleDensity);
%
points = [Shape.vertexPoss, Samples];
Mdl = KDTreeSearcher(points');
numP = size(points, 2);

COLS = knnsearch(Mdl, points', 'K', 37);
COLS = COLS(:, 2:37);
ROWS = kron((1:size(points,2))', ones(1, 36));

ROWS = reshape(ROWS, [1, numP*36]);
COLS = reshape(COLS, [1, numP*36]);

dif = points(:, ROWS) - points(:, COLS);
dis = sqrt(sum(dif.*dif));
tp = sort(dis);
sigma = median(tp(floor(length(tp)/3)));
weights = exp(-dis.*dis/2/sigma/sigma);
G = sparse(double(ROWS), double(COLS), double(weights), numP, numP);
G = max(G, G');
d = full(sum(G));
D = sparse(1:numP, 1:numP, 1./sqrt(d));
L_nor = D*(sparse(1:numP, 1:numP, d) - G)*D;
[u,v] = eigs(L_nor, dimBasis, -1e-10);

numV = size(Shape.vertexPoss, 2);
Basis.eigVecs = u(1:numV, :);
Basis.eigVals = diag(v)';

function [Samples] = mesh_sampling(Shape, sampleDensity)
% add samples within 
numSamples2 = 32768;
numFaces = size(Shape.faceVIds, 2);
p1s = Shape.vertexPoss(:, Shape.faceVIds(1,:));
p2s = Shape.vertexPoss(:, Shape.faceVIds(2,:));
p3s = Shape.vertexPoss(:, Shape.faceVIds(3,:));
e12s = p1s - p2s;
e13s = p1s - p3s;
e23s = p2s - p3s;
as = sum(e12s.*e12s);
bs = sum(e13s.*e13s);
cs = sum(e23s.*e23s);
faceAreas = sqrt(max(0, 2*(as.*bs + as.*cs + bs.*cs) - (as.*as+ bs.*bs + cs.*cs)))/4;
for faceId = 2:numFaces
    faceAreas(faceId) = faceAreas(faceId) + faceAreas(faceId-1);
end

faceAreas = faceAreas/faceAreas(numFaces);
sample_ts = sort(rand(1, numSamples2));
points = zeros(3, numSamples2);

faceId = 1;
for sId = 1 : numSamples2
    while sample_ts(sId) > faceAreas(faceId)
        faceId = faceId + 1;
    end
    p1 = Shape.vertexPoss(:, Shape.faceVIds(1, faceId));
    p2 = Shape.vertexPoss(:, Shape.faceVIds(2, faceId));
    p3 = Shape.vertexPoss(:, Shape.faceVIds(3, faceId));
    r1 = rand(1,1);
    r2 = rand(1,1);
    t1 = (1-sqrt(r1));
    t2 = sqrt(r1)*(1-r2);
    t3 = sqrt(r1)*r2;
    points(:, sId) = t1*p1 + t2*p2 + t3*p3;
end

Mdl = KDTreeSearcher(Shape.vertexPoss');
IDX = knnsearch(Mdl, points');
dis = points - Shape.vertexPoss(:, IDX');
dis = sqrt(sum(dis.*dis));
Samples = points(:, find(dis > sampleDensity));
