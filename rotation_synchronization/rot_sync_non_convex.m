function [R_opt, RR_opt] = rot_sync_non_convex(I, RR, sigma)
% The main function
% I : a 2xN matrix that specifies the input graph,
% the first row indicates indices of the source shape
% and the second row indicates indices of the target shape
%
% RR: a 3x3xN matrix in correspondence with I
%     RR(:,:,i) is a 3x3 orthogonal matrix that specifies the relative
%     roations from object I(1,i) to object I(2,i)
% sigma: a small value. The default value is 1e-2.

% Output: R_opt: 3x3xM matrix encoding the optimized rotation for each
%                object
%         RR_opt: 3x3xN matrix encoding the optimized relative
%         transformation. It is computed from R_opt

if length(sigma) == 0
    sigma = 1e-2;
end

edgeWeights = ones(1, size(I,2));
% Robust rotation synchronization
for iter = 1:30
    [R_opt, lambda] = lap_rot_sync(I, RR, edgeWeights);
    
    if iter < 20
        edgeWeights = weight_maps(I, RR, R_opt, sigma, 1);
    else
        edgeWeights = weight_maps(I, RR, R_opt, sigma, 2);
    end
    fprintf('iteration = %d, lambda = %f.\n', iter, lambda);
end

Tp = reshape(R_opt, [9, size(R_opt,3)]);
RR_opt = helper_matrix_mul(Tp(:,I(2,:)),...
    helper_matrix_transpose(Tp(:,I(1,:))));
RR_opt = reshape(RR_opt, [3,3, size(I,2)]);

function [weights] = weight_maps(I, RR, R_approx, sigma, p)
%
numEdges = size(I, 2);

Ri = reshape(R_approx(:,:, I(1,:)), [9, numEdges]);
Ri = helper_matrix_transpose(Ri);
Rj = reshape(R_approx(:,:, I(2,:)), [9, numEdges]);
Rij = helper_matrix_mul(Rj, Ri);
d = Rij - reshape(RR, [9, numEdges]);
d = sum(d.*d);
sigma2 = sigma^2;
weights = sigma2./(sigma2 + d);

% downweights noisy pairs more when p is big
weights = power(weights, p/2);