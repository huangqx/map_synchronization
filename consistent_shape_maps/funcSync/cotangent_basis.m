function [Basis] = cotangent_basis(mesh, numEigs)
% Compute the Laplacian basis for a triangular mesh
%
numVertices = size(mesh.vertexPoss, 2);
numFaces = size(mesh.faceVIds, 2);

rows = zeros(1, 3*numFaces);
cols = zeros(1, 3*numFaces);
vals = zeros(1, 3*numFaces);

areaWeights = zeros(numVertices, 1);

for faceId = 1 : numFaces
    v1Id = mesh.faceVIds(1, faceId);
    v2Id = mesh.faceVIds(2, faceId);
    v3Id = mesh.faceVIds(3, faceId);
    p1 = mesh.vertexPoss(:, v1Id);
    p2 = mesh.vertexPoss(:, v2Id);
    p3 = mesh.vertexPoss(:, v3Id);
    [cotans, area] = triangle_geometry(p1, p2, p3);
    areaWeights(v1Id) = areaWeights(v1Id) + area/3;
    areaWeights(v2Id) = areaWeights(v2Id) + area/3;
    areaWeights(v3Id) = areaWeights(v3Id) + area/3;
    
    rows(3*faceId - 2) = v2Id;
    cols(3*faceId - 2) = v3Id;
    vals(3*faceId - 2) = cotans(1)/2;
    
    rows(3*faceId - 1) = v1Id;
    cols(3*faceId - 1) = v3Id;
    vals(3*faceId - 1) = cotans(2)/2;
    
    rows(3*faceId) = v1Id;
    cols(3*faceId) = v2Id;
    vals(3*faceId) = cotans(3)/2;
end

L = sparse(rows, cols, vals, numVertices, numVertices);
L = L + L';
L = sparse(1:numVertices, 1:numVertices, sum(L)) - L;
D = sparse(1:numVertices, 1:numVertices, areaWeights);


[eigVecs, eigVals] = eigs(L, D, numEigs, -0.01);
eigVals = diag(eigVals);
[s, ids] = sort(eigVals);

Basis.vertexWeights = areaWeights;
Basis.vecs = eigVecs(:,ids);
Basis.vals = eigVals(ids);


function [cotans, area] = triangle_geometry(p1, p2, p3)
%

e12 = p1 - p2;
l3 = e12'*e12;

e13 = p1 - p3;
l2 = e13'*e13;

e23 = p2 - p3;
l1 = e23'*e23;

area = sqrt(2*(l1*l2 + l2*l3 + l1*l3) - l1*l1 - l2*l2 - l3*l3)/4;


cos1 = min(0.999, max(-0.999, (l2 + l3 - l1)/2/sqrt(l2*l3)));
cos2 = min(0.999, max(-0.999, (l1 + l3 - l2)/2/sqrt(l1*l3)));
cos3 = min(0.999, max(-0.999, (l1 + l2 - l3)/2/sqrt(l1*l2)));

cotans = [cos1/sqrt(1-cos1*cos1);
    cos2/sqrt(1-cos2*cos2);
    cos3/sqrt(1-cos3*cos3)];
