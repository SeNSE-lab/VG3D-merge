function [R,T,varargout] = Get_TransformToStandardOrientation(x0,y0,z0,varargin)
%% function [R,T,{x,y,z}] = Get_TransformToStandardOrientation(x0,y0,z0,{'param_name',param})
% -------------------------------------------------------------------------
% INPUT:
%   (x0,y0,z0) - coordinates of original line object
%       -> If any one coordinate is empty, use a []
%   varargin:
%       'plot' - toggle to plot output
%       'lt' - threshhold for linear region
%       'id_lt' - bypass threshold and specifiy index for linear region
%       'pt' - threshold for planar region
%       'id_pt' - bypass threshold and specify index for planar region
%       'zdir_th' - threshold for determining curvature direction for R(1)
%       'R1_offset' - user-defined offset for zeta, R(1), in degrees
% OUTPUT:
%   R - rotation angles to place in standard orientation
%       [Rx Ry Rz]
%   T - translation to place in standard orientation
%       [Tx Ty Tz]
%   varargout:
%       {x,y,z} - coordinates of line object in standard orientation
% -------------------------------------------------------------------------
% NOTES:
% + Standard orientation is based on definition from RatMap
%       “The Morphology of the Rat Vibrissal Array: A Model 
%       for Quantifying Spatiotemporal Patterns of Whisker-Object Contact”
%       Towal RB*, Quist BW*, Gopal V, Solomon JH, and Hartmann MJZ
%       * authors contributed equally
%       PLoS Computational Biology (2011)    
% + Based on FindOrientationAngles.m from original RatMap code
% -------------------------------------------------------------------------
% HELPER FUNCTIONS:
%   -> createLine3d.m (geom3d toolbox)
%   -> distancePointLine3d.m (geom3d toolbox)
%   -> createPlane.m (geom3d toolbox)
%   -> distancePointPlane (geom3d toolbox)
% -------------------------------------------------------------------------
% Brian Quist
% November 2, 2011

%% Handle inputs
TGL_plot = 0;
R1_offset = 0; % degrees
lt = 0.1; % Linear threshold
pt = 0.1; % Planar threshold
id_lt = NaN;
id_pt = NaN;
zdir_th = 0; % When determing R(1), need to know orientation of curvature
             % This threshold helps set when to flip the curvature
             % which is assumed concave down
if ~isempty(varargin),
   for ii = 1:2:(length(varargin)) 
       switch varargin{ii},
           case 'R1_offset',R1_offset = varargin{ii+1};
           case 'lt',       lt = varargin{ii+1};
           case 'id_lt',    id_lt = varargin{ii+1};
           case 'pt',       pt = varargin{ii+1};
           case 'id_pt',    id_pt = varargin{ii+1};
           case 'plot',     TGL_plot = varargin{ii+1};
           case 'zdir_th',  zdir_th = varargin{ii+1};
           otherwise,
               error('Not a valid input parameter');
       end
   end
end
zdir_th = abs(zdir_th);

%% Handle empty vectors
N = max([length(x0),length(y0),length(z0)]);

% Handle missing data
if isempty(x0), x0 = zeros(1,N); end
if isempty(y0), y0 = zeros(1,N); end
if isempty(z0), z0 = zeros(1,N); end

% Correct vector size
if size(x0,2) == 1, x0 = x0'; end
if size(y0,2) == 1, y0 = y0'; end
if size(z0,2) == 1, z0 = z0'; end

%% Move to (0,0,0)
T = [x0(1) y0(1) z0(1)];
x = x0 - T(1);
y = y0 - T(2);
z = z0 - T(3);

%% Determine linear region
if isnan(id_lt),
    
    maxresidual = zeros(length(x),1);
    for ii = 2:length(x)
        
        referenceline = createLine3d([x(1) y(1) z(1)],[x(ii) y(ii) z(ii)]);
        resline = distancePointLine3d([x(1:ii)',y(1:ii)',z(1:ii)'],referenceline);
        maxresidual(ii) = max(abs(resline));
        
    end
    id_lt = find(maxresidual <= lt,1,'last');
    
end

%% Determine planar region
if isnan(id_pt),
    
    maxresidual = zeros(length(x),1);
    for ii = [2:(id_lt-1),(id_lt+1):length(x)]
        
        refplane = createPlane([x(1) y(1) z(1)],[x(id_lt) y(id_lt) z(id_lt)],...
            [x(ii) y(ii) z(ii)]);
        residuals = zeros(ii,1);
        for jj = 1:ii
            residuals(jj) = distancePointPlane([x(jj) y(jj) z(jj)],refplane);
        end
        maxresidual(ii) = max(abs(residuals));
        
    end
    id_pt = find(maxresidual <= pt,1,'last');
    
end

%% Rz: R(3) ... TH
% Rotate segment to be co-planar with the x-axis
R(3) = atan2(y(id_lt)-y(1),x(id_lt)-x(1))*(180/pi);
if R(3) < 0, R(3) = R(3) + 360; end
[x,y,z] = Get_RotateTranslate(x,y,z,[0 0 -R(3)],[0 0 0]); % -'ve is correct

%% Ry: R(2) ... PHI
% Rotate segment to be aligned with the x-axis 
R(2) = atan2(z(id_lt)-z(1),x(id_lt)-x(1))*(180/pi);
R(2) = -R(2); % Correct negative sign
[x,y,z] = Get_RotateTranslate(x,y,z,[0 -R(2) 0],[0 0 0]); % -'ve is correct

%% Rz: R(1) ... ZETA ... based on planar region
% Determine planar region, then align with xy-plane
pt_normal = planeNormal( createPlane( [x(1) y(1) z(1)], ...
    [x(id_lt) y(id_lt) z(id_lt)],[x(id_pt) y(id_pt) z(id_pt)]));
R(1) = atan2(pt_normal(3),pt_normal(2))*(180/pi);
if isnan(R(1)), R(1) = 0; end

% Adjust R(1) to match RatMap

% normal in first quadrant
if (R(1) >= 0 && R(1) <= 90) && z(end) >= zdir_th
    R(1) = 90 - R(1);
elseif (R(1) >= 0 && R(1) <= 90) && z(end) <= -zdir_th
    R(1) = -(90 + R(1));
    
    % normal in second quadrant
elseif R(1) > 90 && z(end) <= -zdir_th
    R(1) =  -(R(1) - 90);
elseif R(1) > 90 && z(end) >= zdir_th
    R(1) =  270 - R(1);
    
    % normal in third quadrant
elseif R(1) < -90 && z(end) <= -zdir_th
    R(1) = -(270 - abs(R(1)));
elseif R(1) < -90 && z(end) >= zdir_th
    R(1) = abs(R(1)) - 90;
    
    % normal in fourth quadrant
elseif (R(1) < 0 && R(1) >= -90) && z(end) <= -zdir_th
    R(1) = -(90 - abs(R(1)));
elseif (R(1) < 0 && R(1) >= -90) && z(end) >= zdir_th
    R(1) = abs(R(1)) + 90;
end

R(1) = R(1) + R1_offset;
R(1) = -R(1); % Correct w/ -ve sign

% Final Rotation
[x,y,z] = Get_RotateTranslate(x,y,z,[-R(1) 0 0],[0 0 0]); % -'ve is correct

%% Plot Check
if TGL_plot,
    figure;
    plot3(x0,y0,z0,'c.'); hold on;
    plot3(x,y,z,'b.');
    [x1,y1,z1] = Get_RotateTranslate(x,y,z,R,T);
    plot3(x1,y1,z1,'mo');
    grid on;
    xlabel('x'); ylabel('y'); zlabel('z');
    xl = get(gca,'XLim'); yl = get(gca,'YLim'); zl = get(gca,'ZLim');
    plot3([xl(1) xl(2)],[0 0],[0 0],'k-','LineWidth',2);
    plot3([0 0],[yl(1) yl(2)],[0 0],'k-','LineWidth',2);
    plot3([0 0],[0 0],[zl(1) zl(2)],'k-','LineWidth',2);
    % ---
    disp('-----');
    disp(['Rx: ',num2str(R(1))]);
    disp(['Ry: ',num2str(R(2))]);
    disp(['Rz: ',num2str(R(3))]);
    disp(' ');
    disp(['Tx: ',num2str(T(1))]);
    disp(['Ty: ',num2str(T(2))]);
    disp(['Tz: ',num2str(T(3))])    
end

%% Final output
varargout{1} = x;
varargout{2} = y;
varargout{3} = z;