function [G] = fm_shape_graph(Shapes, Para)
% Input arguments:
%       Shapes: input shapes
%       Para.knn: the number of nearest neighbors per object
%       Para.gridDim: the dimension of the grid
% Compute shape descriptors
dess = zeros(Para.gridDim^3, length(Shapes));
for id = 1 : size(dess, 2)
    dess(:, id) = fm_shape_descriptor(Shapes{id}, Para.gridDim);
end
%
similarity = zeros(size(dess,2), size(dess,2));
for id = 1:size(dess, 2)
    for id2 = 1:size(dess, 2)
        similarity(id, id2) = norm(dess(:, id) - dess(:, id2));
    end
end
G = sparse(size(dess, 2), size(dess, 2));
for id = 1 : size(dess, 2)
    [s, neighborIds] = sort(similarity(id,:));
    G(id, neighborIds(2:(Para.knn+1))) = 1;
end
G = max(G, G');