function [x,y] = Get_cesaro2cart(s_total,CC,varargin)
%% function [x,y] = Get_cesaro2cart(s_total,CC,{'param_name',param_value})
% -------------------------------------------------------------------------
% INPUT:
%   s_total - total arc length 
%   CC - Cesaro coefficients such that:
%       k = CC(1)*ds^2 + CC(2)*ds + CC(3)
%   varargin:
%       'N' - number of points
%       'ds' - incremental ds (overrides 's' argument')
% OUTPUT:
%   (x,y) - 2D Cartesian coordinates of resulting line segment
% -------------------------------------------------------------------------
% NOTES:
%   + Based on calculations in cesaro2cart by Joe Solomon
% -------------------------------------------------------------------------
% Brian Quist
% October 5, 2011

%% Constants and variable inputs
N = 50; % 50 points
ds = NaN;

if ~isempty(varargin),
   for ii = 1:2:length(varargin) 
       switch varargin{ii},
           case 'N', N = varargin{ii+1};
           case 'ds', ds = varargin{ii+1}; s = NaN;
           otherwise,
               error('Not a valid input parameter');
       end
   end
end

% Subtract 1 from N to get correct total N:
N = N-1;

%% Setup ds
if isnan(ds),
   ds = ones(1,N).*(s_total/N);    
end
s = [0 cumsum(ds)]; 

%% Compute 2D shape
% Compute change in curvature
k = CC(1).*(s(2:end-1)).^2 + CC(2).*s(2:end-1) + CC(3);
% Angle difference along arc length
dTheta = k.*ds(1:end-1);
theta = cumsum(dTheta);
x = [0 ds(1) ds(1)+cumsum(ds(2:end).*cos(theta))];
y = [0 0 cumsum(ds(2:end).*sin(theta))];