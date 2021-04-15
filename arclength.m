%% * function [s_total,s,ds,x,y] = arclength(x,y,{z,nNodes})
function [s_total,s,ds,x,y] = arclength(x,y,varargin)

%Returns the arc length of curve defined by (x,y,{z}).
% INPUTS:
%   x = x-coordinates
%   y = y-coordinates
%   varargin: (in any order)
%       z - z-coordinates
%       nNodes - number of sample nodes
% OUTPUTS:
%   s_total = total arc length
%   s = arclength at each (x,y,{z}) (starts at zero, ends at s_tot)
%   ds = diff(s)
% Written by Joe Solomon for elasticaPB
% Revised:
%   + March 24, 2009 - Brian Quist - added nNodes sampling capability
%   + March 13, 2010 - Brian Quist - handle (x,y,z) as column
%   + Aug 16, 2012 - Lucie Huet - handle 3D input

% Check varargin for 3D or nNodes defined
z = zeros(size(x));
nNodes = [];
for ii = 1:length(varargin)
    switch length(varargin{ii})
        case length(x), z = varargin{ii};
        case 1,         nNodes = varargin{ii};
    end
end

if size(x,2) == min(size(x)), x = x'; end
if size(y,2) == min(size(y)), y = y'; end
if size(z,2) == min(size(z)), z = z'; end

if ~isempty(nNodes) && nNodes<length(x)
    index = round(1:(length(x)-1)/(nNodes-1):length(x));
    x = x(index);
    y = y(index);
    z = z(index);
end

ds = sqrt((x(2:end)-x(1:end-1)).^2+(y(2:end)-y(1:end-1)).^2+(z(2:end)-z(1:end-1)).^2);
s = [0 cumsum(ds)];
s_total = s(end);