function [x,y,z] = Get_RotateTranslate(x0,y0,z0,R,T)
%% function [x,y,z] = Get_RotateTranslate(x0,y0,z0,R,T)
% -------------------------------------------------------------------------
% INPUT:
%   (x0,y0,z0) - original 3D data points
%       * any of the three inputs can be []
%   R - rotation angles vector (in degrees)
%       [Rx Ry Rz]
%       * can be [] for zero rotation vector
%   T - translation vector
%       [Tx Ty Tz]
%       * can be [] for zero translation vector
% OUTPUT:
%   (x,y,z) - rotated and translated 3D data points
% -------------------------------------------------------------------------
% NOTES:
%   + Based on LOCAL_Calculate3DBasePoints from Get_RatMap
% -------------------------------------------------------------------------
% Brian Quist
% October 6, 2011

%% Structure input data
N = max([length(x0),length(y0),length(z0)]);

% Handle missing data
if isempty(x0), x0 = zeros(1,N); end
if isempty(y0), y0 = zeros(1,N); end
if isempty(z0), z0 = zeros(1,N); end
if isempty(R), R = [0 0 0]; end
if isempty(T), T = [0 0 0]; end

% Correct vector size
if size(x0,2) == 1, x0 = x0'; end
if size(y0,2) == 1, y0 = y0'; end
if size(z0,2) == 1, z0 = z0'; end

%% Setup Rotation matrix

% Determine rotation matrix multiplication factors:
% (Code used to determine the matrix:)
% syms c_x s_x c_y s_y c_z s_z
% B = [1 0 0;0 c_x -s_x; 0 s_x c_x];  % Rotation about x
% C = [c_y 0 s_y; 0 1 0; -s_y 0 c_y]; % Rotation about y
% D = [c_z -s_z 0; s_z c_z 0; 0 0 1]; % Rotation about z
% A = D*C*B;
% A =
% [ c_y*c_z, c_z*s_x*s_y - c_x*s_z, s_x*s_z + c_x*c_z*s_y]
% [ c_y*s_z, c_x*c_z + s_x*s_y*s_z, c_x*s_y*s_z - c_z*s_x]
% [    -s_y,               c_y*s_x,               c_x*c_y]

%% Rotate and Translate

d2r = pi/180;

% Matrix multiplication
c_x = cos(R(1)*d2r);
s_x = sin(R(1)*d2r);
c_y = cos(R(2)*d2r);
s_y = sin(R(2)*d2r);
c_z = cos(R(3)*d2r);
s_z = sin(R(3)*d2r);
A = [...
    c_y*c_z, c_z*s_x*s_y - c_x*s_z, s_x*s_z + c_x*c_z*s_y; ...
    c_y*s_z, c_x*c_z + s_x*s_y*s_z, c_x*s_y*s_z - c_z*s_x; ...
    -s_y,    c_y*s_x,               c_x*c_y];

% Rotate
x = x0.*A(1,1) + y0.*A(1,2) + z0.*A(1,3);
y = x0.*A(2,1) + y0.*A(2,2) + z0.*A(2,3);
z = x0.*A(3,1) + y0.*A(3,2) + z0.*A(3,3);

% Translate
x = x + T(1);
y = y + T(2);
z = z + T(3);