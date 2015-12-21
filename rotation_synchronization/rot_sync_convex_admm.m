function [R_opt, RR_opt] = rot_sync_convex_admm(I, RR)
% The main function
% I : a 2xN matrix that specifies the input graph,
% the first row indicates indices of the source shape
% and the second row indicates indices of the target shape
%
% RR: a 3x3xN matrix in correspondence with I
%     RR(:,:,i) is a 3x3 orthogonal matrix that specifies the relative
%     roations from object I(1,i) to object I(2,i)
% 
% Output: R_opt: 3x3xM matrix encoding the optimized rotation for each
%                object
%         RR_opt: 3x3xN matrix encoding the optimized relative
%         transformation. It is computed from R_opt



