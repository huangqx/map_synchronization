function [Basis] = fm_man_made_basis(Shape, dimBasis)
% Compute the functional basis for a given 3D model
numFaces = size(Shape.faceVIds, 2);
%
G_face = compute_face_graph(Shape);
G_face2 = G_face*G_face;
%
[faceCenters, faceNors] = compute_face_normals(Shape);
%
[rows, cols, vals] = find(G_face2);
ids = find(rows < cols);
rows = rows(ids)';
cols = cols(ids)';

% Compute the weights between adjacent faces
edgeLen = faceCenters(:, rows) - faceCenters(:, cols);
edgeLen = sqrt(sum(edgeLen.*edgeLen));
sigmaLen = median(edgeLen);
%weightLen = exp(-(edgeLen.*edgeLen)/2/sigmaLen/sigmaLen);
weightLen = exp(-edgeLen/sigmaLen);
%
nor1 = faceNors(:, rows);
nor2 = faceNors(:, cols);
edgeAngle = acos(min(1, abs(sum(nor1.*nor2))));
sigmaAngle = max(0.2, median(edgeAngle));
weightAngle = exp(-edgeAngle/sigmaAngle);
weight = weightLen.*weightAngle;
%
Adj_face = sparse(rows, cols, weight, numFaces, numFaces);
Adj_face = Adj_face + Adj_face';
%
d = sum(Adj_face);
L_face = sparse(1:length(d), 1:length(d), d) - Adj_face;
[U,V] = eigs(L_face, dimBasis, 1e-16);
%
[s, ids] = sort(diag(V)');
Basis.eigVals = diag(V(ids,ids))';
Basis.eigVecs = U(:, ids);

function [flag] = is_neighbor(vids1, vids2)
flag = length(intersect(vids1, vids2)) == 2;

function [G_face] = compute_face_graph(Shape)
%
numV = size(Shape.vertexPoss, 2);
numF = size(Shape.faceVIds, 2);
%
ADJ_vf = sparse(Shape.faceVIds, ones(3,1)*(1:numF), ones(3, numF),...
    numV, numF);
G_face = sparse(numF, numF);
for vId = 1 : numV
    faceIds = find(ADJ_vf(vId, :));
    for i = 1:length(faceIds)
        fid1 = faceIds(i);
        for j = (i+1):length(faceIds)
            fid2 = faceIds(j);
            if is_neighbor(Shape.faceVIds(:, fid1), Shape.faceVIds(:, fid2)) == 1
                G_face(fid1, fid2) = 1;
            end
        end
    end
end
G_face = max(G_face, G_face');
%
h = 10;

function [faceCenters, faceNors] = compute_face_normals(Shape)
%
Pos1 = Shape.vertexPoss(:, Shape.faceVIds(1, :));
Pos2 = Shape.vertexPoss(:, Shape.faceVIds(2, :));
Pos3 = Shape.vertexPoss(:, Shape.faceVIds(3, :));
%
faceCenters = (Pos1 + Pos2 + Pos3)/3;
e12 = Pos1 - Pos2;
e13 = Pos1 - Pos3;
%
faceNors = cross(e12, e13);
norms = sqrt(sum(faceNors.*faceNors));
faceNors = faceNors./(ones(3,1)*norms);

