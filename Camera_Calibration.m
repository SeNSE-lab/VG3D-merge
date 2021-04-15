%% Camera_Calibration
% -------------------------------------------------------------------------
% NOTES:
% + Requires: 
%   Camera Calibration Toolbox for Matlab
%   http://www.vision.caltech.edu/bouguetj/calib_doc/
% -------------------------------------------------------------------------
% Brian Quist
% November 10, 2011
global CAMERA_TOP CAMERA_FRONT

%% [A] Setup paths:
%CAMERA_FRONT = 'F:\20111010_Calibration_PROC\Camera No.1_C001H001S0001'; %#ok<*NASGU>
%CAMERA_TOP = 'F:\20111010_Calibration_PROC\Camera No.2_C002H001S0001';

%% [B] Run calibration for FRONT camera (Camera #1)
close all;
cd(CAMERA_FRONT);

% Steps:
data_calib;     % (1) 'Image names'
click_calib;    % (2) 'Extract grid corners'
go_calib_optim; % (3) 'Calibration'
saving_calib;   % (4) 'Save'

% ... OR can access via: calib_gui_normal;

%% [C] Run calibration for TOP camera (Camera #2)
close all;
cd(CAMERA_TOP);

% Steps:
data_calib;     % (1) 'Image names'
click_calib;    % (2) 'Extract grid corners'
go_calib_optim; % (3) 'Calibration'
saving_calib;   % (4) 'Save'

% ... OR can access via: calib_gui_normal;

%% [D] Move data to higher directory
% FRONT Camera
clear; global CAMERA_FRONT %#ok<*REDEF>
cd(CAMERA_FRONT);
load Calib_Results;
cd ..; clear CAMERA_FRONT
save Calib_Results_FRONT

% TOP Camera
clear; global CAMERA_TOP
cd(CAMERA_TOP);
load Calib_Results;
cd ..; clear CAMERA_TOP
save Calib_Results_TOP

%% [E] Run Stereo calibration
% * Copy calibration files to parent folder of FRONT/TOP images
% * Add _FRONT and _TOP to data filenames
% * LEFT = FRONT
% * RIGHT = TOP

% Steps:
load_stereo_calib_files;    % (1) Load left and right calibration files
go_calib_stereo;            % (2) Run stereo calibration
saving_stereo_calib;        % (3) Save stereo calib results

% ... OR can access via: stereo_gui;

%% [F] Compact Stereo calibration
Camera_Calibration_CompactStereo;