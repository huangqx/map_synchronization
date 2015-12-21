function [Rt] = helper_matrix_transpose(R)
% This function computes the transpose of rotation matrices in the
% vectorized form

Rt = R;
Rt(2,:) = R(4,:);
Rt(3,:) = R(7,:);
Rt(4,:) = R(2,:);
Rt(6,:) = R(8,:);
Rt(7,:) = R(3,:);
Rt(8,:) = R(6,:);