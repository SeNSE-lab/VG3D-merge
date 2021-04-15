function [x,y,z,summary_PT] = Merge3D(xf,yf,xt,yt,count,varargin)
% disp('3D_Merging.m')
global summary_PT
%% [x,y,z] = Merge3D(xf,yf,xt,yt,{'param_name',param})
% -------------------------------------------------------------------------
% INPUT:
%   (xf,yf) - (x,y) points corresponding to line in image A ... FRONT (Left) Camera
%   (xt,yt) - (x,y) points corresponding to line in image B ... TOP (Right) Camera
%   varargin
%       'plot' - toggle plotting
%       'camera_path' - path to camera calibration parameters
%           -> Ouput of Camera_Calibration_CompactStereo.m 
%       'wm_opts' - cell array of 3D worm options
%       'bp_opts' - cell array of base point options
% -------------------------------------------------------------------------
% NOTES:
% + Uses stereo_triangulation from:
%   Camera Calibration Toolbox for Matlab
%   http://www.vision.caltech.edu/bouguetj/calib_doc/
%
% + TOP and FRONT camera views should have same object orientation
%   e.g. top-right of checkboard is top-right in both views
% -------------------------------------------------------------------------
% Brian Quist
% November 7, 2011

%% Handle inputs 
TGL_plot = 1;
camera_path = [];
wm_opts = {};
bp_opts = {};
if ~isempty(varargin),
   for ii = 1:2:(length(varargin)) 
       switch varargin{ii},
           case 'wm_opts', wm_opts = varargin{ii+1};
           case 'bp_opts', bp_opts = varargin{ii+1};
           case 'camera_path', camera_path = varargin{ii+1};
           case 'plot', TGL_plot = varargin{ii+1};
           otherwise,
               error('Not a valid input parameter');
       end
   end
end

%% Setup Camera calibration parameters
if isempty(camera_path)
    aaaaaa = which('Calib_Results_stereo.mat');
    disp('Loading camera calibration from:');
    disp(aaaaaa);
    load Calib_Results_stereo
    clear aaaaaa
else
    aaaaaa = pwd;
    cd(camera_path);
    load Calib_Results_stereo
    cd(aaaaaa);
    clear aaaaaa
end    

%% Compute best-fitting base-point
[bp_x,bp_y,bp_z,er,Front_projBP,Top_projBP] = Fit_3dBasepoint(xf,yf,xt,yt,...
    'A_proj',{fc_left,cc_left,kc_left,alpha_c_left}, ...
    'B_proj',{fc_right,cc_right,kc_right,alpha_c_right}, ...
    'A2B_transform',{om,T}, ...
    'Plot_Final',TGL_plot,bp_opts{:});

%%  JAE addition 140925
%   Aims to fix basepoint mismatches on the spot
%   You need to be sure that the xf,yf vectors correspond to the same view
%   as REF.bp_Ap, and xt,yt to REF_bp_Bp
BP_er_thresh_low = 2;
BP_er_thresh_high = 100;
comments = {'We re-fitting that basepoint for you, right now!','Your basepoint is a little off, let me get that for you...','Basepoints are for the birds...and rats...and yours is off.  Fixing!','BP FIX MODE!!!','BP FIX, FTW!!'};
if er > BP_er_thresh_low && er < BP_er_thresh_high
    disp(comments{randsample(1:length(comments),1)})
   [xf,yf,xt,yt] = shall_not_be_denied(xf,yf,xt,yt,Front_projBP,Top_projBP);
   [bp_x,bp_y,bp_z] = Fit_3dBasepoint(xf,yf,xt,yt,...
    'A_proj',{fc_left,cc_left,kc_left,alpha_c_left}, ...
    'B_proj',{fc_right,cc_right,kc_right,alpha_c_right}, ...
    'A2B_transform',{om,T}, ...
    'Plot_Final',TGL_plot,bp_opts{:});
else disp(sprintf('er=%d',er))
end

%% Use a set basepoint from attempts
% chosen_point = 2;
% load BP_reffy
% disp(['new BP = ',num2str(REF.attempted_3d_bp{end}')])
% bp_x = REF.attempted_3d_bp{chosen_point}(1);
% bp_y = REF.attempted_3d_bp{chosen_point}(2);
% bp_z = REF.attempted_3d_bp{chosen_point}(3);

%% Fit 3D Worm
%% Ellis control clause
% if abs(length(xf)-length(xt)) < 200
[x,y,z,PT] = Fit_3dWorm(xf,yf,xt,yt, ...
    'BP',[bp_x,bp_y,bp_z], ...
    'A_proj',{fc_left,cc_left,kc_left,alpha_c_left}, ...
    'B_proj',{fc_right,cc_right,kc_right,alpha_c_right}, ...
    'A2B_transform',{om,T}, ...
    'Plot_Final',TGL_plot,wm_opts{:});
% else
%     PT=[];
%     x=0;
%     y=0;
%     z=0;
% end

close all

summary_PT{count} = PT;
