## Merge3D

This repository includes code used in Bush et al. 2021 *in prep* to create a 3D whisker shape from whiskers tracked in two separate camera views. 

The primary entry point to the merging process is the file `mergeLoop.m`. This function takes the following inputs:

- tw - The top tracked whisker struct
- fw - The front tracked whisker struct
- calibInfo - a [1x10] cell array of calibration info (for more info see https://docs.opencv.org/master/d9/db7/tutorial_py_table_of_contents_calib3d.html): {fc_left,cc_left,kc_left,alpha_c_left,fc_right,cc_right,kc_right,alpha_c_right,om,T}
- tracked_3D_fileName - the output filename

and outputs a struct with the 3D points(x,y,z) of the whisker for each recorded frame.

`tw` and `fw` are N-length structs with fields x and y, where N is the number of video frames from which the points were tracked. These structs follow directly from the output of "Whisk": https://www.janelia.org/open-science/whisk-whisker-tracking

The quality of the merge depends strongly on accurate estimation of the basepoint, as the merge considers the first point for each frame of `tw` and `fw` as a starting point. The merge then  guesses a point in 3D space that is the next point along the whisker. It then reprojects that 3D point into the 2D camera views and calculates the error between this reprojected 2D point and the tracked 2D whisker. The merge repeats this reprojection and attempts to minimize the error of the reprojected point. Once this error has been minimized, the next point along the whisker is estimated from the previously found point, until the entire whisker is tracked.

Accurate camera calibration using a checkerboard of known dimension is crucial for accurate 3D reconstruction

Preprocessing the 2D whisker shapes to be spatially smooth is also advised.