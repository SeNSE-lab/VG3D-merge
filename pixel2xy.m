function [x,y] = pixel2xy(bw_image,varargin)
%% function [x,y] = pixel2xy(bw_image,{'setting_name',setting})
% -------------------------------------------------------------------------
% INPUT:
%   bw_image - binary image to convert to (x,y) coordinates
%   varargin:
%       'convert_type' - either 1 {default} or 0 to be converted to (x,y)
%       'coordinate_frame' - either 'pixel' {default} or 'cart'
%       'plot' - toggle output plot, 0 {default} or 1
% OUTPUT:
%   (x,y) - coordinates of extracted image
% -------------------------------------------------------------------------
% NOTES:
%   .
% -------------------------------------------------------------------------
% Brian Quist
% October 19, 2011

%% Setup inputs
TGL_PltFinal = 0; 
TGL_ConvertType = 1;
TGL_CoordinateFrame = 'pixel';
if ~isempty(varargin),
   for ii = 1:2:(length(varargin)) 
       switch varargin{ii},
           case 'convert_type', TGL_ConvertType = varargin{ii+1};
           case 'coordinate_frame', TGL_CoordinateFrame = varargin{ii+1};
           case 'plot', TGL_PltFinal = varargin{ii+1};
           otherwise,
               error('Not a valid input parameter');
       end
   end
end

%% Setup outputs
if TGL_ConvertType == 1,
    n_pts = sum(sum(bw_image));
else
    n_pts = sum(sum(~bw_image));
end
x = zeros(n_pts,1);
y = zeros(n_pts,1);

%% Convert bw_image to (x,y) coordinates
ctr = 1;
for ii = 1:size(bw_image,2)
    iy = find(bw_image(:,ii) == TGL_ConvertType);
    if ~isempty(iy),
        for z = 1:length(iy)
            x(ctr) = ii; 
            switch TGL_CoordinateFrame,
                case 'pixel', y(ctr) = iy(z); 
                case 'cart',  y(ctr) = size(bw_image,1) - iy(z); 
                otherwise
                    error('Bad coordinate frame selection');
            end
            ctr = ctr+1;
        end
    end
end

%% Plot
if TGL_PltFinal,
    figure;
    imshow(bw_image); hold on;
    plot(x,y,'b.');
end