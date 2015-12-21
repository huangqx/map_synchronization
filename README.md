This repository (will) provide implementations of several map synchronization algorithms that compute maps among a collection of shapes. The code is released under the MIT license and can be used for any purpose with proper attribution. The code accompanies the following paper, which should be cited in publications that use the provided modules:


*****************************************************************
Folder 'io' provides matlab codes for saving and loading shapes in wavefront obj format *****************************************************************

*****************************************************************
Folder “manmade_align” implements the joint alignment method described in

"Fine-Grained Semi-Supervised Labeling of Large Shape Collections", Qixing Huang, Hao Su, and Leonidas Guibas. SIGGRAPH ASIA' 13.

It takes a collection of shapes of the same class as input, and aligns them in a common space. 

The shapes are assumed to have given upright orientations. Two example datasets (100 chairs and 100 cars) are provided in folder data/

The code proceeds in two steps. The first step optimizes an affine transformation to align all the input shapes. The second optimizes a free-from deformation to align the input shapes. 

% Step I
Shapes_aff = man_made_joint_affine_main(Shapes_in, Para_align);

% Step II
Shapes_ffd = man_made_joint_ffd_main(Shapes_aff, Para_align);


*****************************************************************

*****************************************************************
Folder “consistent_shape_maps” implements the joint alignment method described in

"Consistent Shape Maps via Semidefinite Programming."Qixing Huang and Leonidas Guibas. Computer Graphics Forum, Volume 32, Issue 5, Proc. Eurographics Symposium on Geometry Processing (SGP), 2013. (Best Paper Award).

It takes a collection of shapes of the same class and noisy pair-wise maps as input, and output consistent maps among them.

% load data set
load('data\csm_human.mat');
% Alternatively, you can go to 'consistent_shape_maps\pointSync'
% and call [Data] = load_dataset(foldername) to load a dataset in 
% 'data\csm'

% load parameters
load ('data\para_csm.mat');

% run the program
[opt_maps] = csm_main_func(Data, Para_csm);


*****************************************************************


Currently, this repository is still being constantly updated. I will provide more implementations soon. 

