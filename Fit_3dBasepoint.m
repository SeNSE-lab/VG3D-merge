function [x,y,z,varargout] = Fit_3dBasepoint(AAx,AAy,BBx,BBy,varargin)
% disp(sprintf('\nFit3dBasepoint.m\n'))
%% function [x,y,z,P] = Fit_3dBasepoint(AAx,AAy,BBx,BBy,{'setting_name',setting})
% -------------------------------------------------------------------------
% INPUT:
%   (AAx,AAy) - (x,y) points corresponding to line in image A ... FRONT (Left) Camera
%   (BBx,BBy) - (x,y) points corresponding to line in image B ... TOP (Right) Camera
%   
%   varargin:
%       
%   * Required Settings:
%       'A_proj' - Camera settings for image A
%           -> Formatted as: {fc,cc,kc,alpha_c}
%       'B_proj' - Camera settings for image B
%           -> Formatted as: {fc,cc,kc,alpha_c}
%       'A2B_transform' - Transformation matrix components for A->B
%           -> Formatted as: {om,T}
%
%   * Optional settings:
%       'Plot_Final' - 1/0 to show intermediate fitting
%       'N' - Number of points to use in base-point approximation
%
% OUTPUT:
%   (x,y,z) - best-fit base-point in 3D coordinates
%
%   varargout:
%       P - Projections
%               .bp_A - selected base-points in A to fit
%               .bp_Ap - projection of best-fit base-point in view A
%               .bp_B - selected base-points in B to fit
%               .bp_Bp - projection of best-fit base-point in view B
% -------------------------------------------------------------------------
% NOTES:
%   + Assumes first points of (AAx,AAy) and (BBx,BBy) are the base
%   + TOP and FRONT camera views should have same object orientation
%     e.g. top-right of checkboard is top-right in both views
% -------------------------------------------------------------------------
% Brian Quist 
% October 26, 2011
% Revised:
%   12/12/2011 - Matt Graff correction to uf/ut swap for correct output
global REF

%% Handle inputs 
TGL_PltFinal = 1;
N = 1;
A_camera = {ones(2,1),zeros(2,1),zeros(5,1),0};
B_camera = {ones(2,1),zeros(2,1),zeros(5,1),0};
A2B_transform = {zeros(3,1),zeros(3,1)};
if ~isempty(varargin),
    for ii = 1:2:(length(varargin))
        switch varargin{ii},
            case 'Plot_Final', TGL_PltFinal = varargin{ii+1};
            case 'N', N = varargin{ii+1};
            case 'A_proj', A_camera = varargin{ii+1};
            case 'B_proj', B_camera = varargin{ii+1};
            case 'A2B_transform', A2B_transform = varargin{ii+1};
            otherwise,
                error('Not a valid input parameter');
        end
    end
end

% Check size if inputs
if size(AAx,2) == 1, AAx = AAx'; end
if size(AAy,2) == 1, AAy = AAy'; end
if size(BBx,2) == 1, BBx = BBx'; end
if size(BBy,2) == 1, BBy = BBy'; end

% Setup global reference parameter
REF.bp_A = [AAx(1:N);AAy(1:N)];
REF.bp_B = [BBx(1:N);BBy(1:N)];
REF.A_camera = A_camera;
REF.B_camera = B_camera;
REF.A2B_transform = A2B_transform;
REF.bp_Ap = [NaN;NaN];
REF.bp_Bp = [NaN;NaN];

% JAE addition 140330
REF.attempted_bp_A = {};
REF.attempted_bp_B = {};
REF.attempted_3d_bp = {};

%% Guess 3D point location
[XL,XR] = stereo_triangulation ...
    (REF.bp_A(:,end),REF.bp_B(:,end),A2B_transform{1},A2B_transform{2}, ...
    A_camera{1},A_camera{2},A_camera{3},A_camera{4}, ...
    B_camera{1},B_camera{2},B_camera{3},B_camera{4}); %#ok<NASGU>

%% Search for best-fit
options = optimset('tolx',1e-6); % fminsearch options
[q_final,er] = fminsearch(@LOCAL_FindBestPoint,[XL(1) XL(2) XL(3)],options); % fit by minimizing error w.r.t Cartesian coords.
% [q_right_final,er] = fminsearch(@LOCAL_FindBestPoint_Jimmy_right,[XL(1) XL(2) XL(3)],options); % fit by minimizing error w.r.t Cartesian coords.
% [q_left_final,er] = fminsearch(@LOCAL_FindBestPoint_Jimmy_left,[q_right_final(1) q_right_final(2) q_right_final(3)],options);

%% Output final
% This uses the front view to set the x coordinate
x = q_final(1);
y = q_final(2);
z = q_final(3);

P.bp_A = REF.bp_A;
P.bp_Ap = REF.bp_Ap;
P.bp_B = REF.bp_B;
P.bp_Bp = REF.bp_Bp;
% varargout{1} = P;

%%  JAE addition 140925
%   For in-line basepoint fix
varargout{1} = er;
varargout{2} = REF.bp_Ap;
varargout{3} = REF.bp_Bp;

%% Final Plot
if TGL_PltFinal,
    figure;
    set(gcf,'Position',[20 250 800 400])
    set(gcf,'Name','Basepoint Calibration Results','NumberTitle','off');
    % ---
    subplot(1,2,1);
    plot(AAx,AAy,'k.'); hold on;
    plot(REF.bp_A(1,:),REF.bp_A(2,:),'bo');
    plot(REF.bp_Ap(1),REF.bp_Ap(2),'r*');
    legend({'Projection';'Fitted base point(s)';'Final fitted point'}, ...
        'Location','SouthWest','FontSize',8);
    axis equal;
    title('FRONT (Left) View');
    % ---
    subplot(1,2,2);
    plot(BBx,BBy,'k.'); hold on;
    plot(REF.bp_B(1,:),REF.bp_B(2,:),'bo');
    plot(REF.bp_Bp(1),REF.bp_Bp(2),'r*');
    axis equal;
    title('TOP (Right) View');
    % ---
    text(0.05,0.99,{ ...
        ['x (XL): ',num2str(q_final(1))]; ...
        ['y (XL): ',num2str(q_final(2))]; ...
        ['z (XL): ',num2str(q_final(3))]; ...
        ['------ ']; ...
        ['er    : ',num2str(er)]}, ...
        'FontSize',8,'Units','normalized', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','top'); %#ok<NBRAK>
    drawnow;

end


function e = LOCAL_FindBestPoint_Jimmy_right(q)
%% function e = LOCAL_FindBestPoint(q)
global REF

% Compute Camera A projection
[uf,vf] = Get_3DtoCameraProjection(q(1),q(2),q(3), ...
    'proj',REF.A_camera);

% Convert 3D point in A coordinate frame to B coordinate frame
% ->  Y = rigid_motion(X,om,T)
r = rigid_motion(q',REF.A2B_transform{1},REF.A2B_transform{2});

% Compute Camera B projection
[ut,vt] = Get_3DtoCameraProjection(r(1),r(2),r(3), ...
    'proj',REF.B_camera);

% REF output
REF.bp_Ap = [uf;vf];
REF.bp_Bp = [ut;vt];

% Compute error
er_f = sum(sqrt((REF.bp_A(1,:)-uf).^2 + (REF.bp_A(2,:)-vf).^2));
er_t = sum(sqrt((REF.bp_B(1,:)-ut).^2 + (REF.bp_B(2,:)-vt).^2));
e = er_t;

function e = LOCAL_FindBestPoint_Jimmy_left(q)
%% function e = LOCAL_FindBestPoint(q)
global REF

% Compute Camera A projection
[uf,vf] = Get_3DtoCameraProjection(q(1),q(2),q(3), ...
    'proj',REF.A_camera);

% Convert 3D point in A coordinate frame to B coordinate frame
% ->  Y = rigid_motion(X,om,T)
r = rigid_motion(q',REF.A2B_transform{1},REF.A2B_transform{2});

% Compute Camera B projection
[ut,vt] = Get_3DtoCameraProjection(r(1),r(2),r(3), ...
    'proj',REF.B_camera);

% REF output
REF.bp_Ap = [uf;vf];
REF.bp_Bp = [ut;vt];

% JAE addition 140330
REF.attempted_bp_A{length(REF.attempted_bp_A)+1} = [uf;vf];
REF.attempted_bp_B{length(REF.attempted_bp_B)+1} = [ut;vt];
% % Guess 3D point location AGAIN -- JAE addition 140331
% [nXL,nXR] = stereo_triangulation ...
%     ([uf;vf],[ut;vt],REF.A2B_transform{1},REF.A2B_transform{2}, ...
%     REF.A_camera{1},REF.A_camera{2},REF.A_camera{3},REF.A_camera{4}, ...
%     REF.B_camera{1},REF.B_camera{2},REF.B_camera{3},REF.B_camera{4}); %#ok<NASGU>
% REF.attempted_3d_bp{length(REF.attempted_3d_bp)+1} = nXL;

% Compute error
er_f = sum(sqrt((REF.bp_A(1,:)-uf).^2 + (REF.bp_A(2,:)-vf).^2));
er_t = sum(sqrt((REF.bp_B(1,:)-ut).^2 + (REF.bp_B(2,:)-vt).^2));
e = er_t;


function e = LOCAL_FindBestPoint(q)
%% function e = LOCAL_FindBestPoint(q)
global REF

% Compute Camera A projection
[uf,vf] = Get_3DtoCameraProjection(q(1),q(2),q(3), ...
    'proj',REF.A_camera);

% Convert 3D point in A coordinate frame to B coordinate frame
% ->  Y = rigid_motion(X,om,T)
r = rigid_motion(q',REF.A2B_transform{1},REF.A2B_transform{2});

% Compute Camera B projection
[ut,vt] = Get_3DtoCameraProjection(r(1),r(2),r(3), ...
    'proj',REF.B_camera);

% REF output
REF.bp_Ap = [uf;vf];
REF.bp_Bp = [ut;vt];

% JAE addition 140330
REF.attempted_bp_A{length(REF.attempted_bp_A)+1} = [uf;vf];
REF.attempted_bp_B{length(REF.attempted_bp_B)+1} = [ut;vt];

% Guess 3D point location AGAIN -- JAE addition 140331
[nXL,nXR] = stereo_triangulation ...
    ([uf;vf],[ut;vt],REF.A2B_transform{1},REF.A2B_transform{2}, ...
    REF.A_camera{1},REF.A_camera{2},REF.A_camera{3},REF.A_camera{4}, ...
    REF.B_camera{1},REF.B_camera{2},REF.B_camera{3},REF.B_camera{4}); %#ok<NASGU>
REF.attempted_3d_bp{length(REF.attempted_3d_bp)+1} = nXL;

% Compute error
er_f = sum(sqrt((REF.bp_A(1,:)-uf).^2 + (REF.bp_A(2,:)-vf).^2));
er_t = sum(sqrt((REF.bp_B(1,:)-ut).^2 + (REF.bp_B(2,:)-vt).^2));
e = er_f + er_t;

