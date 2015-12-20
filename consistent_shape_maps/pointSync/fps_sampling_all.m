function [SAMPLE] = fps_sampling_all(Data, numSamples)
% Perform sampling of all the shapes
for shapeId = 1:length(Data.shapes)
    SAMPLE{shapeId} = fps_sampling(Data.shapes{shapeId}, numSamples);
    fprintf('Finished sampling shape %d\n', shapeId);
end

function [SAMPLE] = fps_sampling(shape, numSamples)
% Perform furthest point sampling to generate uniform samples on a given
% model

% SAMPLE will be used in the algorithm

% Generate the graph used to compute shortest path
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

% vertex indices of samples
SAMPLE.sampleIds = zeros(1, numSamples); 
% the closest sample index of each vertex
SAMPLE.closestSampleIds = zeros(1, numPoints);
% The matrix that stores pair-wise distances between samples and mesh
% points
SAMPLE.distMat = zeros(numSamples, numPoints);

% Compute the nextSampleId
dist = graphshortestpath(G, 1);
[maxDis, nextSampleId] = max(dist);
% distance vector for picking the next sample
distVec = ones(1, numPoints)*1000;

% Iterative sampling
for i = 1:numSamples
    if i <= length(shape.featurePointIds)
%        nextSampleId = shape.featurePointIds(i);
    end
    SAMPLE.sampleIds(i) = nextSampleId;
    dist = graphshortestpath(G, nextSampleId);
    SAMPLE.distMat(i, :) = dist;
    ids = find(dist < distVec);
    distVec(ids) = dist(ids);
    SAMPLE.closestSampleIds(ids) = i;
    [maxDis, nextSampleId] = max(distVec);
end
fprintf('density = %f.\n', maxDis);

% Compute the distance matrix between feature points and mesh vertices
% feaDisMat = zeros(length(shape.featurePointIds), numPoints);
% for i = 1:length(shape.featurePointIds)
%     dist = graphshortestpath(G, shape.featurePointIds(i));
%     feaDisMat(i,:) = dist;
% end
