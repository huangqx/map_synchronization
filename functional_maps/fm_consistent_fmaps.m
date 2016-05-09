function [fmaps] = fm_consistent_fmaps(Data, Para)
% Input arguments:
%       Data.shapes:             the input shapes      
%       Data.basis:              the precomputed basis
%       Para.lambda_regularize:  the regularization prior for pair-wise
%                                functional maps
%       Para.lambda_consistency: the consistency we enforce on the shape
%                                network
%
