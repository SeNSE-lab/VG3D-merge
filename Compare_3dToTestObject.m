function [e_al,e,stats] = Compare_3dToTestObject(x,y,z,obj,varargin)
%% function [e_al,e,stats] = Compare_3dToTestObject(x,y,z,obj,{'param_name',param})
% -------------------------------------------------------------------------
% INPUT:
%   (x,y,z) - 3D reconstructed coordinates of reference object
%   obj - reference object type
%           -> 'n','s1','s2','w'
%   varargin
%       'opts_TSO_3d' - Options for tranform to standard orientation for 3D data
%       'opts_TSO_2d' - Options for tranform to standard orientation for 2D data
%       'plot' - toggle plotting
%       'text' - toggle text
% OUTPUT:
%   e_al - error difference in arc lengths
%   e - error w.r.t. 3D object
%   stats - struct of error fit stats
%       .max - maximum error
%       .min - minimum error
%       .median - median error
%       .std - standard deviation of error
%       .N - number of points in calculation
% -------------------------------------------------------------------------
% NOTES:
%   + (setting) and its required 2D file:
%       ('n') : CLIPPED_20111010_Nitinol.tif
%       ('s1'): CLIPPED_20111010_SSWire_Type1.tif
%       ('s2'): CLIPPED_20111010_SSWire_Type2_rotated.tif
%       ('w') : CLIPPED_20111010_Whisker.tif
%   + These files should be in the path already before running the function
%     e.g. addpath('C:\ ... Image Location Folder');
% -------------------------------------------------------------------------
% Brian Quist
% November 7, 2011

%% [0] Handle variable inputs
TGL_plot = 1;
TGL_text = 1;
opts_TSO_3d = {};
opts_TSO_2d = {}; %#ok<NASGU>
if ~isempty(varargin),
   for ii = 1:2:(length(varargin)) 
       switch varargin{ii},
           case 'opts_TSO_3d', opts_TSO_3d = varargin{ii+1};
           case 'opts_TSO_2d', opts_TSO_2d = varargin{ii+1}; %#ok<NASGU>
           case 'plot', TGL_plot = varargin{ii+1};
           case 'text', TGL_text = varargin{ii+1};
           otherwise,
               error('Not a valid input parameter');
       end
   end
end

%% [2D] Pixels2mm
% Scans taken at 1200 dpi
% 1 / ((1200 pixels/inch) * (1 inch/25.4 mm))
p2mm = 1/(1200 * (1/25.4)); % millimeters / pixel

%% [2D] Load image
switch obj
    case 'n', 
        q0 = imread('CLIPPED_20111010_Nitinol.tif');
        N = 0.01;
        RL = 0.01;
        opts_TSO_2d = {'R1_offset',90};
    case 's1',    
        q0 = imread('CLIPPED_20111010_SSWire_Type1.tif');
        N = 0.001;
        RL = 0.05;
        opts_TSO_2d = {'R1_offset',180};
    case 's2',    
        q0 = imread('CLIPPED_20111010_SSWire_Type2_rotated.tif');
        N = 0.001;
        RL = 0.05;
        opts_TSO_2d = {'R1_offset',180};
    case 'w', 
        q0 = imread('CLIPPED_20111010_Whisker.tif');
        N = 0.03;
        RL = 0.01;
        opts_TSO_2d = {'R1_offset',180};
    otherwise
        error('Wrong object type');
end

%% [2D] Extract
qt = im_cleanborder(q0,'fill_v',256,'bwidth',10);
bwt = im_bw(qt,'invert',true,'bw_thresh',0.95,'max_prop','Area');
[xt_2d,yt_2d] = pixel2xy(bwt);
xt_2d = xt_2d.*p2mm;
yt_2d = yt_2d.*p2mm;

%% [2D] Find centerline and arc length
[xc_2d,yc_2d] = Get_centerline(xt_2d,yt_2d,'N',round(length(xt_2d)*N), ...
    'R', sqrt((xt_2d(end)-xt_2d(1))^2+(yt_2d(end)-yt_2d(1))^2)*RL);
[xc_2d,yc_2d] = equidist(xc_2d,yc_2d);

%% [2D] Place in standard orientation (based on RatMap)
[~,~,x_2d,y_2d,z_2d] = Get_TranformToStandardOrientation(xc_2d,yc_2d,[],'plot',0,opts_TSO_2d{:});

%% [3D] Place in standard orientation (based on RatMap)
[~,~,x_3d,y_3d,z_3d] = Get_TranformToStandardOrientation(x,y,z,'plot',0,opts_TSO_3d{:});

%% [COMPARE] 3D shapes
if TGL_plot,
   figure; 
   set(gcf,'Name','2D reference object and 3D object','NumberTitle','off');
   % ---
   plot3(x_2d,y_2d,z_2d,'b.'); hold on;
   plot3(x_3d,y_3d,z_3d,'co');
   legend({'Reference object';'3D object'},'Location','SouthWest','FontSize',8);
   xlabel('x'); ylabel('y'); zlabel('z');
   xl = get(gca,'XLim'); yl = get(gca,'YLim'); zl = get(gca,'ZLim');
   plot3([xl(1) xl(2)],[0 0],[0 0],'k-','LineWidth',2);
   plot3([0 0],[yl(1) yl(2)],[0 0],'k-','LineWidth',2);
   plot3([0 0],[0 0],[zl(1) zl(2)],'k-','LineWidth',2);
   grid on;
   axis equal;    
end

%% [COMPARE] Arc length
s_2d = arclength3d(x_2d,y_2d,z_2d);
s_3d = arclength3d(x_3d,y_3d,z_3d);
e_al = abs(s_2d-s_3d);
if TGL_text,
    disp('=============');
    disp('+ Arclength:');
    disp([' 2D: ',num2str(s_2d)]);
    disp([' 3D: ',num2str(s_3d)]);
    disp([' ----------']); %#ok<NBRAK>
    disp(['     ',num2str(e_al)]);
end

%% [COMPARE] Find maximum error
[e,stats] = Get_ShapeError(x_2d,y_2d,z_2d,x_3d,y_3d,z_3d,'plot',TGL_plot);
if TGL_text,
    disp('=============');
    disp('Error stats:');
    disp('-------------');
    disp('  Absolute:');
    disp([' Max: ',num2str(stats.max)]);
    disp([' mean&std: ',num2str(stats.mean), ...
        ' +/- ',num2str(stats.std)]);
    disp('-------------');
    disp('  % (100) of s (2D):');
    disp([' Max: ',num2str(stats.max/stats.s_total*100)]);
    disp([' mean&std: ',num2str(stats.mean/stats.s_total*100), ...
        ' +/- ',num2str(stats.std/stats.s_total*100)]);
end
