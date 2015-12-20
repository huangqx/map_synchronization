function [feaDisMat] = eval_precompt(Data)

for i = 1:length(Data.shapes)
    feaDisMat{i} = compt_fea_dis(Data.shapes{i});
end

function [feaDisMat] = compt_fea_dis(shape)

numPoints = size(shape.vertexPoss, 2);
sIds = [shape.faceVIds(1,:), shape.faceVIds(2,:), shape.faceVIds(3,:)];
tIds = [shape.faceVIds(2,:), shape.faceVIds(3,:), shape.faceVIds(1,:)];
G = sparse(double(sIds), double(tIds), ones(1, length(sIds)), numPoints, numPoints);
G = G + G';

if 1
    G = G*G; % Use a 2-ring neighborhood graph
end

[sIds, tIds, vals] = find(G);
difVec = shape.vertexPoss(:, sIds) - shape.vertexPoss(:, tIds);
difVec = sqrt(sum(difVec.*difVec));
G = sparse(sIds, tIds, double(difVec), numPoints, numPoints);

feaDisMat = zeros(length(shape.featurePointIds), numPoints);
for i = 1:length(shape.featurePointIds)
    dist = graphshortestpath(G, shape.featurePointIds(i));
    feaDisMat(i,:) = dist;
end