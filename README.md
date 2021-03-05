# ConnScript_RAT

Scripts for preprocessing resting-state fMRI data in Don Clusterio.

Other sites/computers may use it, just check the definitions for the atlas and how to call fsl ($FSLDIR), ants and afni.

The scripts change the voxel-size to avoid sub-milimeter dimensions (currently x10, the atlas is also changed accordingly).
Note that some lines are aimed at reorienting the images as the Atlas.
Also some variables are defined at the beginning of the file, check accordingly.
Although we like BIDS for the folder/files naming, we are not fully bids-compatible.

This is just a simple script.
We provide no warranty that it will work with your data, it may be a good starting point.

Enjoy and comment,

Sarael Alcauter

GitHub version, march 5, 2021 
