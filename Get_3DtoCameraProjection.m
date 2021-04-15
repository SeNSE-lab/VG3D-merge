function [u,v] = Get_3DtoCameraProjection(x,y,z,varargin)
%% function [u,v] = Get_3DtoCameraProjection(x,y,z,{'param_name',param})
% -------------------------------------------------------------------------
% INPUT:
%   (x,y,z) - Object in 3D coordinates (row vectors)
%   varargin
%       'proj' - projection plane
%           -> 'XY' - xy-plane ... 'Aerial' view
%           -> 'YZ' - yz-plane ... 'Front' view
%           -> 'ZX' - zx-plane ... 'Side' view
%           -> camera plane ... input cell of camera parameters
%               + Input structure from: {fc,cc,kc,alpha_c}
% OUTPUT:
%   (u,v) - 2D projection
% -------------------------------------------------------------------------
% NOTES:
% + Camera parameters {fc,cc,kc,alpha_c} use output from:
%       Camera Calibration Toolbox for Matlab
%       http://www.vision.caltech.edu/bouguetj/calib_doc/
% + For math of camera projection see:
%       http://opencv.willowgarage.com/documentation/
%       camera_calibration_and_3d_reconstruction.html
% -------------------------------------------------------------------------
% Brian Quist
% October 26, 2011

%% Handle inputs
proj = 'XY';
c_params = [];
if ~isempty(varargin),
   for ii = 1:2:(length(varargin)) 
       switch varargin{ii},
           case 'proj', proj = varargin{ii+1};
           otherwise,
               error('Not a valid input parameter');
       end
   end
end

%% Check structure of (x,y,z)
if size(x,2) == 1, x = x'; end
if size(y,2) == 1, y = y'; end
if size(z,2) == 1, z = z'; end

%% Compute projection
if ~ischar(proj), 
    c_params = proj;
    proj = 'C';
end
switch proj
    case 'XY', u = x; v = y; % Aerial View
    case 'YZ', u = y; v = z; % Front View
    case 'ZX', u = x; v = z; % Side View
    case 'C', 
        % Camera View
        if ~isempty(c_params)
            
            % Setup camera parameters
            fc = c_params{1};
            cc = c_params{2};
            kc = c_params{3};
            alpha_c = c_params{4};
            
            % Assume camera is at origin of global frame
            om = [0;0;0];
            T = [0;0;0];
            
            % Compute projection
            xp = project_points2([x;y;z],om,T,fc,cc,kc,alpha_c);
            u = xp(1,:);
            v = xp(2,:);      
            
        else
            error('Camera parameters not supplied');
        end        
    otherwise,
        error('Bad choice of projection view');
end