%% Settings_VG3D_20130313_T09_F0001F4800
% -------------------------------------------------------------------------
% Settings file for elastica3D processing
% -------------------------------------------------------------------------
% DATA SPECIFICS:
%   + Data author: Chris Schroeder
%   + Date taken: 3/13/2013
%   + Trial: 9
%   + Frames: 1-4800
% -------------------------------------------------------------------------
% Lucie Huet
% April 24, 2013

%% 00. Path names and labeling
% ** CHANGE WITH EACH DATASET **

% Frame sequence
PT.Frames = [1 4800];

% Save path location:
PT.save = 'C:\WHISKER\MATLAB\BendData';

% Data path location:
% PT.dataname = 'Camera No.2_C002H001S0009';
PT.data = 'C:\WHISKER\MATLAB\BendData'; % TOP

% Data TAG:
PT.TAG = ['_120313_T09_F',num2str(PT.Frames(1),'%04.0f'),'F',num2str(PT.Frames(2),'%04.0f')];

%% Run_Extract_2D_WSKR
% --> Run: Setup_PEG_Top_v2.m
% =========================================================================
% ** CHANGE WITH EACH DATASET: **

% /////////////////////////////////////////////////////////////////////////
% ** Change as needed (should not vary much across datasets) **
% PT.ImageOrientation = 'S'; % Direction that nose points

%% Run_Process_3D_WSKR
% ** CHANGE WITH EACH DATASET: **

% /////////////////////////////////////////////////////////////////////////
% ** Change as needed (should not vary much across datasets) **
PT.REF = NaN; % put specific reference frames here
% REF can either be individual frame numbers which will be applied to the
% entire contact they are a part of, OR it can be an array the length of C
% defining values at each frane so you can specify different reference
% frames within a single contact

PT.pix2m = 0.000058; % conversion to change pixels to meters

PT.pct_planar = 0.70; % what percentage of the undeflected whisker we want 
% planar (for calculating ZETA)
% Knutsen 2008 says 0.70 planar
% Towal 2011 says loosely 0.65 or 0.50 strictly

% also check sign of TH_world!!!


%% Elastica3D
PT.E3D_E = 3.3e9; % Young's Modulus, 3.3 GPa
PT.E3D_rbase = 100e-6; % radius of base, 100 micrometers
PT.E3D_taper = 15; % taper ratio (base to tip)
PT.filt = Inf; % a number from 0 to 1 for cutoff frequency, 1 is half the sample rate
% for no filter, put PT.filt to Inf

