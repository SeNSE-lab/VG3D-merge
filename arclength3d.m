function [s_total,s,ds,x,y,z] = arclength3d(x,y,z,varargin)
%% function [s_total,s,ds,x,y,z] = arclength3d(x,y,z,{nNodes})
% -------------------------------------------------------------------------
% + Returns the arc length of curve defined by (x,y,z)
% -------------------------------------------------------------------------
% INPUT:
%   x = x-coordinates
%   y = y-coordinates
%   z = z-coordinates
%   varargin:
%       {1}: nNodes - number of sample nodes
% OUTPUT:
%   s_total = total arc length
%   s = arclength at each (x,y,z) (starts at zero, ends at s_tot)
%   ds = diff(s)
% -------------------------------------------------------------------------
% NOTES:
%   + Based originally on arclength.m by Joe Solomon
% -------------------------------------------------------------------------
% Brian Quist
% October 31, 2011

%% Handle x,y,z as a column
if size(x,2) == min(size(x)), x = x'; end
if size(y,2) == min(size(y)), y = y'; end
if size(z,2) == min(size(z)), z = z'; end

%% Redo indexes
if ~isempty(varargin),
    nNodes = varargin{1};
    if nNodes < length(x)
        index = round(1:(length(x)-1)/(nNodes-1):length(x));
        x = x(index);
        y = y(index);
        z = z(index);
    end
end

%% Compute ds and arclength
ds = sqrt( ...
    (x(2:end)-x(1:end-1)).^2 + ... 
    (y(2:end)-y(1:end-1)).^2 + ...
    (z(2:end)-z(1:end-1)).^2);
s = [0 cumsum(ds)];
s_total = s(end);