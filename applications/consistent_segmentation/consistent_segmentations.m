function [Segs] = consistent_segmentations(Data, Para)
% Input arguments:
%       Data.shapes:   input shapes
%       Data.basis:    the basis per shape
%       Data.fmaps:    the pre-computed pair-wise fmaps
% Output argument:
%       Segs:          the face segment ids



function [Unary] = unary_potential(Shape, Basis, Para)
% Encoding the image segmentation cues using a quadratic form

