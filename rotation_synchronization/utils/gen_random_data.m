function [I, RR] = gen_random_data(n, p_density, q_error)
% Generate a random graph with density p_density and the probability of
% wrong rotations is q_error
F = rand(n,n) < p_density;
[rows, cols, vals] = find(F);

ids = find(rows < cols);
I = [rows(ids)'; cols(ids)'];

numEdges = size(I, 2);

RR = zeros(3,3, numEdges);

flags = rand(1,numEdges) < q_error;

for eId = 1 : numEdges
    if flags(eId) == 1
        RR(:,:,eId) = rs_rand3rot();
    else
        RR(:,:,eId) = eye(3);
    end
end
