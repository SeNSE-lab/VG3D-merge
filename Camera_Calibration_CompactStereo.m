%% Camera_Calibration_CompactStereo
% -------------------------------------------------------------------------
% NOTES:
% + Start with stereo calibration from:
%   Camera Calibration Toolbox for Matlab
%   http://www.vision.caltech.edu/bouguetj/calib_doc/
%
% + With completed calibration, will save critical values
% -------------------------------------------------------------------------
% Brian Quist
% October 21, 2011

%% Clear out
clear; close all; clc;

%% Load stero file
load Calib_Results_stereo

%% Save stero file critical components
save Calib_Results_stereo_compact ...
    om T ...
    fc_left  cc_left  kc_left  alpha_c_left  ...
    fc_right cc_right kc_right alpha_c_right