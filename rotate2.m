%% * function [x2,y2] = rotate2(x1,y1,theta(radians),origin)
function [x2,y2] = rotate2(x1,y1,theta,origin)

%Rotates points (x1,y1) by angle theta about specified origin.
% INPUTS:
%   x1 = input x-coordinates
%   y1 = input y-coordinates
%   theta = angle to rotate by (counter-clockwise)
%   origin = origin about which to rotate (vector of length 2)
% OUTPUTS:
%   x2 = output x-coordinates
%   y2 = output y-coordinates
% -------------------------------------------------------------------------
% Originally by: Joe Solomon
%   Revised: 2010/04/27 - Brian Quist - Add variable length (x,y,theta)
%            2010/05/28 - Brian Quist - Perform vector operation to improve speed
%            2010/06/03 - Brian Quist - Fixed small bug for vector operation

if nargin == 3
    origin = [0 0];
end

if size(x1,1) == 1, x1 = x1'; end
if size(y1,1) == 1, y1 = y1'; end

x1 = x1 - origin(1);
y1 = y1 - origin(2);

if length(theta) == 1,
    x2 = (x1.*cos(theta) + y1.*(-sin(theta))) + origin(1);
    y2 = (x1.*sin(theta) + y1.*cos(theta)) + origin(2);
    x2 = x2';
    y2 = y2';
else
    q = zeros(length(x1),2);
    for ii = 1:length(x1)
        q(ii,:) = [x1(ii) y1(ii)]*[cos(theta(ii)) sin(theta(ii)); -sin(theta(ii)) cos(theta(ii))];
    end
    x2 = q(:,1)' + origin(1);
    y2 = q(:,2)' + origin(2);
end
