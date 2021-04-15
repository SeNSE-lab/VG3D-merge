function [CC,S,R] = Fit_CameraViewsTo3D(AAx,AAy,BBx,BBy,varargin)
disp('Fit_CameraViewsTo3D.m')
%% function [CC,S,R] = Fit_CameraViewsTo3D(AAx,AAy,BBx,BBy,{'setting_name',setting})
% -------------------------------------------------------------------------
% INPUT:
%   (AAx,AAy) - (x,y) points corresponding to line in image A ... FRONT (Left) Camera
%   (BBx,BBy) - (x,y) points corresponding to line in image B ... TOP (Right) Camera
%   varargin:
%       'Plot_Steps' - 1/0 to show intermediate fitting
%       'CC_0' - Guess for Cesaro coefficients
%       'S_0'  - Guess for arc length
%       'R_0'  - Guess for rotation matrix
%       'A_proj' - projection selection for image A (for Get_3DtoCameraProjection)
%       'B_proj' - projection selection for image B (for Get_3DtoCameraProjection)
% OUTPUT:
%   CC - Cesaro coefficients such that:
%       k = CC(1)*ds^2 + CC(2)*ds + CC(3)
%   S - Arc length of 3D segment
%   R - rotation angles vector (in degrees), 
%       R: [Rx Ry Rz]
% -------------------------------------------------------------------------
% NOTES:
%   + Assumes first point of (AAx,AAy) and (BBx,BBy) are the same in 3D space
% -------------------------------------------------------------------------
% Brian Quist 
% October 6, 2011
global S_min TGL_PltSteps AAX AAY BBX BBY Axc Ayc Bxc Byc A_proj B_proj

%% Handle inputs
TGL_PltSteps = 1; 
TGL_PltFinal = 1;
CC_0 = NaN;
S_0 = NaN;
R_0 = NaN;
N = 50;     % Number of points in 3D segment (see also Get_cesaro2cart.m)
A_proj = 'YZ';
B_proj = 'XY';
if ~isempty(varargin),
   for ii = 1:2:(length(varargin)) 
       switch varargin{ii},
           case 'Plot_Steps', TGL_PltSteps = varargin{ii+1};
           case 'Plot_Final', TGL_PltFinal = varargin{ii+1};
           case 'CC_0', CC_0 = varargin{ii+1};
           case 'S_0', S_0 = varargin{ii+1};
           case 'R_0', R_0 = varargin{ii+1};
           case 'N', N = varargin{ii+1};
           case 'A_proj', A_proj = varargin{ii+1};
           case 'B_proj', B_proj = varargin{ii+1};
           otherwise,
               error('Not a valid input parameter');
       end
   end
end

%% Initial Guess: CC
if isnan(CC_0),
    CC_0 = [0 0 0];
end

%% Initial Guess: R
if isnan(R_0),
    Rz = atan2(AAy(2)-AAy(1),AAx(2)-AAx(1))*(180/pi);
    Ry = -atan2(BBy(2)-BBy(1),BBx(2)-BBx(1))*(180/pi);
    Rx = 0;
    R_0 = [Rx Ry Rz];
end

%% Initial Guess: S
if isnan(S_0),
    % Compute arc length from each view
    sA = sum(sqrt((AAx(2:end)-AAx(1:end-1)).^2+((AAy(2:end)-AAy(1:end-1)).^2)));
    sB = sum(sqrt((BBx(2:end)-BBx(1:end-1)).^2+((BBy(2:end)-BBy(1:end-1)).^2)));
    % Guess is maximum arc length found
    S_0 = max([sA,sB]);
end
S_min = S_0;

%% Compute reference matricies
AAX = repmat(AAx',1,N);
AAY = repmat(AAy',1,N);
BBX = repmat(BBx',1,N);
BBY = repmat(BBy',1,N);

%% Search for best-fit
options = optimset('tolx',1e-6); % fminsearch options
q_final = fminsearch(@LOCAL_FindBestCurve,[CC_0 R_0 S_0],options); % fit by minimizing error w.r.t Cartesian coords.

%% Output final
CC = [q_final(1) q_final(2) q_final(3)];
R = [q_final(4) q_final(5) q_final(6)];
S = q_final(7);

%% Final Plot
if TGL_PltFinal,
    LOCAL_PlotResults(AAX,AAY,Axc,Ayc,BBX,BBY,Bxc,Byc,CC,S,R);
end

function e = LOCAL_FindBestCurve(q)
%% function e = LOCAL_FindBestCurve(q)
global S_min TGL_PltSteps AAX AAY BBX BBY Axc Ayc Bxc Byc A_proj B_proj

% Setup Inputs:
CCs = [q(1) q(2) q(3)];
R = [q(4) q(5) q(6)];
S = q(7);

% Constrain S minimum length
% if S < S_min, S = S_min; end

% Compute 3D shape
[x,y] = Get_cesaro2cart(S,CCs);
[x,y,z] = Get_RotateTranslate(x,y,[],R,[]);

% Compute reference projections
[Axc,Ayc] = Get_3DtoCameraProjection(x,y,z,'proj',A_proj);
[Bxc,Byc] = Get_3DtoCameraProjection(x,y,z,'proj',B_proj);

% Compute error for view A: Front (Left)
AXs = repmat(Axc,size(AAX,1),1);
AYs = repmat(Ayc,size(AAX,1),1);
E = sqrt((AXs-AAX).^2 + (AYs-AAY).^2);
eA = min(E);

% Compute error for view B: Top (Right)
BXs = repmat(Bxc,size(BBX,1),1);
BYs = repmat(Byc,size(BBX,1),1);
E = sqrt((BXs-BBX).^2 + (BYs-BBY).^2);
eB = min(E);

% Sum errors to get final error
e = sum(eA + eB) + 1/S;

% Plot
if TGL_PltSteps,
    LOCAL_PlotResults(AAX,AAY,Axc,Ayc,BBX,BBY,Bxc,Byc,CCs,S,R);
end

function LOCAL_PlotResults(AAX,AAY,Axc,Ayc,BBX,BBY,Bxc,Byc,CCs,S,R)
figure(101); clf(101);
set(101,'Position',[20 250 800 400])
% ---
subplot(1,2,1);
plot(AAX(:,1),AAY(:,1),'ko-'); hold on;
plot(Axc,Ayc,'r.-');
axis equal;
title('FRONT (Left) View');
% ---
subplot(1,2,2);
plot(BBX(:,1),BBY(:,1),'ko-'); hold on;
plot(Bxc,Byc,'r.-');
axis equal;
title('TOP (Right) View');
% ---
text(0.05,0.99,{ ... 
    ['CC(1): ',num2str(CCs(1))]; ...
    ['CC(2): ',num2str(CCs(2))]; ...
    ['CC(3): ',num2str(CCs(3))]; ...
    ['R(1) : ',num2str(R(1))]; ...
    ['R(2) : ',num2str(R(2))]; ...
    ['R(3) : ',num2str(R(3))]; ...
    ['S    : ',num2str(S)]}, ...
    'FontSize',8,'Units','normalized', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top');
drawnow;
