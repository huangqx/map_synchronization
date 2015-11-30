This repository (will) provide implementations of several map synchronization algorithms that compute maps among a collection of shapes.

*****************************************************************
folder manmade_align/ implements the joint alignment method described in

"Fine-Grained Semi-Supervised Labeling of Large Shape Collections", Qixing Huang, Hao Su, and Leonidas Guibas. SIGGRAPH ASIA' 13.

It takes a collection of shapes of the same class as input, and aligns them in a common space. 

The shapes are assumed to have given upright orientations. Two example datasets (100 chairs and 100 cars) are provided in folder data/

The code proceeds in two steps. The first step optimizes an affine transformation to align all the input shapes. The second optimizes a free-from deforamtion to align the input shapes. 

% Step I
Shapes_aff = man_made_joint_affine_main(Shapes_in, Para_align);
% Step II
Shapes_ffd = man_made_joint_ffd_main(Shapes_aff, Para_align);


*****************************************************************

Currently, this repository is still being constantly updated. I will provide more implementations soon. 

