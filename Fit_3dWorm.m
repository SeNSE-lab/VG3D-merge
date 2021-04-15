function [wx,wy,wz,PT] = Fit_3dWorm(AAX,AAY,BBX,BBY,varargin)
% disp(sprintf('\nFit_3dWorm.m\n'))
%% function [wx,wy,wz] = Fit_3dWorm(AAX,AAY,BBX,BBY,{'setting_name',setting})
% -------------------------------------------------------------------------
% INPUT:
%   (AAX,AAY) - (x,y) points corresponding to line in image A
%               ... FRONT (Left) Camera (e.g. [y,z] projection)
%   (BBX,BBY) - (x,y) points corresponding to line in image B
%               ... TOP (Right) Camera (e.g. [x,y] projection)
%   varargin:
%
%       + Fit Related:
%       'BP' - (x,y,z) 3D point of base. If not provided, assumes (0,0,0)
%       'DS' - incremental segment length to fit 3D object
%       'N' - number of sub-divisions for DS to include when fitting
%       'ER_thresh' - Error (in pixels) to limit for each projection
%       'PreviousFit' - {wx wy wz} from previous frame (to improve speed)
%                       * note: input is a cell vector
%                       * note: pass in {[]} for first frame if calling this setting
%       'PreviousFit_GuessOnly' - [1]: Use previous fit as guess only
%                                 [0]: Use previous fit if within tolerance {default}
%
%       + fminsearch Related:
%       'set_fminsearch' - {'param',param_setting} settings for optimset
%
%       + Camera Related:
%       'A_proj' - projection selection for image A (for Get_3DtoCameraProjection)
%           + either 'YZ' or cell array of camera parameters {fc,cc,kc,alpha_c}
%       'B_proj' - projection selection for image B (for Get_3DtoCameraProjection)
%           + either 'XY' or cell array of camera parameters {fc,cc,kc,alpha_c}
%       'A2B_transform' - Transformation matrix components for A->B
%           + Formatted as: {om,T}
%           + Required if A_proj and B_proj are camera projections
%
%       + Toggles
%       'Plot_Steps' - 1/0 to show intermediate fitting
%       'Plot_Final' - 1/0 to show final fit
%       'spline' - 1/0 to toggle spline smoothing of final 3D data
%
% OUTPUT:
%   (wx,wy,wz) - 3D points fit by the 3D worm
% -------------------------------------------------------------------------
% NOTES:
%   + A2B_transform is required if A_proj and B_proj are camera projections
%   + TOP and FRONT camera views should have same object orientation
%     e.g. top-right of checkboard is top-right in both views
% -------------------------------------------------------------------------
% Brian Quist
% March 2, 2012
global TGL_PltSteps DS PT mAAX mAAY mBBX mBBY A_proj B_proj A2B

DS = []; PT = []; mAAX = []; mAAY = []; mBBX = []; mBBY = []; A2B = [];%#ok<NASGU>

%% Handle inputs
TGL_PltSteps = 0;
TGL_PltFinal = 0;
TGL_spline = 0;
TGL_previousfit = 0;
TGL_preivousfit_guessonly = 0;

% Defaults
% -------------------------------------------------------------------------
DS = 1;
PT.N = 5; % Number of subdivisions for DS to test % Default was 5
ER_thresh = 5; % Pixels error in projections to accept
BP = [0 0 0];
A_proj = 'YZ';
B_proj = 'XY';
A2B_transform = {zeros(3,1),zeros(3,1)};

% Previous fit guess
previousfit = [];
Hwx = NaN;
Hwy = NaN;
Hwz = NaN;

% fminsearch settings
fmin_settings = {'maxfunevals',500};
exitflag = NaN;
exitval = NaN;

% User inputs
% -------------------------------------------------------------------------
if ~isempty(varargin),
    for ii = 1:2:(length(varargin))
        switch varargin{ii},
            case 'spline',  TGL_spline = varargin{ii+1};
            case 'DS',      DS = varargin{ii+1};
            case 'BP',      BP = varargin{ii+1};
            case 'N',       PT.N = varargin{ii+1};
            case 'W',       PT.W = varargin{ii+1};
            case 'ER_thresh', ER_thresh = varargin{ii+1};
            case 'A_proj',  A_proj = varargin{ii+1};
            case 'B_proj',  B_proj = varargin{ii+1};
            case 'A2B_transform', A2B_transform = varargin{ii+1};
            case 'Plot_Steps', TGL_PltSteps = varargin{ii+1};
            case 'Plot_Final', TGL_PltFinal = varargin{ii+1};
            case 'PreviousFit', previousfit = varargin{ii+1};
            case 'set_fminsearch', fmin_settings = varargin{ii+1};
            case 'PreviousFit_GuessOnly', TGL_preivousfit_guessonly = varargin{ii+1};
            otherwise,
                error('Not a valid input parameter');
        end
    end
end
A2B = A2B_transform; % Make the transform global

% Transfer previousfit
if ~isempty(previousfit),
    Hwx = previousfit{1};
    Hwy = previousfit{2};
    Hwz = previousfit{3};
    TGL_previousfit = true;
    clear previousfit
end

%% Project base-point into each view
% Update BP if provided (overwrites 'BP' setting)
if TGL_previousfit,
    Hii = 1; % Index of history point to use
    BP(1) = Hwx(Hii);
    BP(2) = Hwy(Hii);
    BP(3) = Hwz(Hii);
    % Increment index for next fit
    Hii = Hii+1;
end

% Compute reference projections: A_proj
[BP_Ax,BP_Ay] = Get_3DtoCameraProjection(BP(1),BP(2),BP(3),'proj',A_proj);

% Convert (x,y,z) from A_frame to B_frame
% ->  Y = rigid_motion(X,om,T)
r = rigid_motion([BP(1);BP(2);BP(3)],A2B{1},A2B{2});

% Compute reference projections; B_proj
[BP_Bx,BP_By] = Get_3DtoCameraProjection(r(1,:),r(2,:),r(3,:),'proj',B_proj);

% Reference output
PT.BP_A = [BP_Ax;BP_Ay];
PT.BP_B = [BP_Bx;BP_By];

clear r BP_Ax BP_Ay BP_Bx BP_By

%% Search 3D space
flag = true;
flip = true;

% Setup base points
PT.wx = BP(1);
PT.wy = BP(2);
PT.wz = BP(3);
PT.BP = [BP(1) BP(2) BP(3)];

% Setup Axc
if TGL_previousfit,
    PT.Axc = PT.BP_A(1);
    PT.Ayc = PT.BP_A(2);
    PT.Bxc = PT.BP_B(1);
    PT.Byc = PT.BP_B(2);
end

% Check size
if size(AAX,1) == 1, AAX = AAX'; end
if size(AAY,1) == 1, AAY = AAY'; end
if size(BBX,1) == 1, BBX = BBX'; end
if size(BBY,1) == 1, BBY = BBY'; end

% Exit check vectors
CHK_A = zeros(size(AAX));
CHK_B = zeros(size(BBX));
clip = NaN;

while flag,
    
    % Short-circuit search if previous fit data provided
    % ---------------------------------------------------------------------
    if TGL_previousfit,
        
        % Compute 3D projection
        % -----------------------------------------------------------------
        % Compute reference projections: A_proj
        [Axc,Ayc] = Get_3DtoCameraProjection( ...
            Hwx(Hii),Hwy(Hii),Hwz(Hii),'proj',A_proj);
        
        % Convert (x,y,z) from A_frame to B_frame
        % ->  Y = rigid_motion(X,om,T)
        r = rigid_motion([Hwx(Hii);Hwy(Hii);Hwz(Hii)],A2B{1},A2B{2});
        
        % Compute reference projections; B_proj
        [Bxc,Byc] = Get_3DtoCameraProjection(r(1,:),r(2,:),r(3,:),'proj',B_proj);
        
        PT.Axc = [PT.Axc; Axc];
        PT.Ayc = [PT.Ayc; Ayc];
        PT.Bxc = [PT.Bxc; Bxc];
        PT.Byc = [PT.Byc; Byc];
        
        % Compute Error
        % -----------------------------------------------------------------
        % Stop if error from single projection is too large
        eA = sqrt( (Axc-AAX).^2 + (Ayc-AAY).^2 );
        eB = sqrt( (Bxc-BBX).^2 + (Byc-BBY).^2 );
        if min(eA) > ER_thresh,
            disp(['eA error (',num2str(min(eA)),') is too large'])
            TGL_previousfit = false;
        elseif min(eB) > ER_thresh,
            disp(['eB error (',num2str(min(eB)),') is too large'])
            TGL_previousfit = false;
        end
        
        % Append data
        % -----------------------------------------------------------------
        if TGL_previousfit,
            PT.wx = [PT.wx Hwx(Hii)];
            PT.wy = [PT.wy Hwy(Hii)];
            PT.wz = [PT.wz Hwz(Hii)];
            % ---
            PT.BP_A = [Axc; Ayc];
            PT.BP_B = [Bxc; Byc];
            % ---
            Hii = Hii+1;
            if Hii > length(Hwx), flag = false; end
        else
            % Start full Fit_3dWorm function
        end
        
        % Toggle for PreviousFit_GuessOnly
        if TGL_preivousfit_guessonly,
            % Stop using previous data and use only as a guess
            TGL_previousfit = false;
        end
    end
    
    % If previousfit data is within the threshold, do not run fit algorithm
    % ---------------------------------------------------------------------
    if ~TGL_previousfit,
        
        % Setup search radius
        % -----------------------------------------------------------------
        % Radius guess
        switch flip,
            case 1,
                PT.RR = max([ ...
                    sqrt( (max(AAX)-min(AAX)).^2 + (max(AAY)-min(AAY)).^2 ), ...
                    sqrt( (max(BBX)-min(BBX)).^2 + (max(BBY)-min(BBY)).^2 )])/10;
            case 2,
                RA = sqrt( ...
                    (PT.Axc(end)-PT.Axc(end-PT.N)).^2 + ...
                    (PT.Ayc(end)-PT.Ayc(end-PT.N)).^2);
                RB = sqrt( ...
                    (PT.Bxc(end)-PT.Bxc(end-PT.N)).^2 + ...
                    (PT.Byc(end)-PT.Byc(end-PT.N)).^2);
                PT.RR = max([RA RB]);
            otherwise
        end
        
        % Reference data
        
        % JAE kludge addition 140401 - work around basepoint for B (top)
        % camera being coerced by increasing PT.RR thresh
        % PT.RR = 100;
        
        % -----------------------------------------------------------------
        % Restrict reference data to circular band
        DA = sqrt((AAX-PT.BP_A(1)).^2 + (AAY-PT.BP_A(2)).^2);
        goodA = logical(DA <= PT.RR  );
        DB = sqrt((BBX-PT.BP_B(1)).^2 + (BBY-PT.BP_B(2)).^2);
        goodB = logical(DB <= PT.RR );
        
        % Setup reference data
        mAAX = repmat(AAX(goodA)',PT.N,1);
        mAAY = repmat(AAY(goodA)',PT.N,1);
        mBBX = repmat(BBX(goodB)',PT.N,1);
        mBBY = repmat(BBY(goodB)',PT.N,1);
        
        % Note data checked
        CHK_A(goodA) = 1;
        CHK_B(goodB) = 1;
        
        if ~TGL_preivousfit_guessonly,
            % Setup guess
            % --------------------------------------------------------------
            switch flip
                case 1
                    DA = sqrt( ...
                        (AAX(goodA)-PT.BP_A(1)).^2 + ...
                        (AAY(goodA)-PT.BP_A(2)).^2);
                    id_A = find(max(DA) == DA,1,'first');
                    tAx = AAX(goodA); tAy = AAY(goodA);
                    vA = [tAx(id_A);tAy(id_A)];
                    clear DA id_A tAx tAy
                    
                    DB = sqrt( ...
                        (BBX(goodB)-PT.BP_B(1)).^2 + ...
                        (BBY(goodB)-PT.BP_B(2)).^2);
                    id_B = find(max(DB) == DB,1,'first');
                    tBx = BBX(goodB); tBy = BBY(goodB);
                    vB = [tBx(id_B);tBy(id_B)];
                    
                    clear DB id_B tBx tBy
                otherwise
                    % View A vector guess
                    vA = LOCAL_ClockGuess(AAX(goodA),AAY(goodA), ...
                        PT.Axc(end-PT.N:end),PT.Ayc(end-PT.N:end),PT.RR);
                    % View B vector guess
                    vB = LOCAL_ClockGuess(BBX(goodB),BBY(goodB), ...
                        PT.Bxc(end-PT.N:end),PT.Byc(end-PT.N:end),PT.RR);
            end
            
            % Stereo-triangulation to guess next 3D point location
            ZZ = stereo_triangulation( ...
                vA,vB, ...
                A2B_transform{1},A2B_transform{2}, ...
                A_proj{1},A_proj{2},A_proj{3},A_proj{4}, ...
                B_proj{1},B_proj{2},B_proj{3},B_proj{4});
        else
            % USE PREVIOUS GUESS
            % --------------------------------------------------------------
            
            % Compute 3D projection
            % -----------------------------------------------------------------
            ZZ(1) = Hwx(Hii);
            ZZ(2) = Hwy(Hii);
            ZZ(3) = Hwz(Hii);
            
            % Increment index for next fit
            Hii = Hii+1;
            
            % Check Hii length
            if Hii > length(Hwx),
                TGL_preivousfit_guessonly = false;
            end
        end
        
        % Compute guess for q0
        R = sqrt( (ZZ(2)-PT.BP(2))^2 + (ZZ(1)-PT.BP(1))^2 );
        % (1) Estimate PHI angle (-'ve b.c. of right hand rule)
        q0(1) = - atan2(ZZ(3)-PT.BP(3),R)*(180/pi);
        % (2) Estimate THETA angle
        q0(2) = atan2(ZZ(2)-PT.BP(2), ZZ(1)-PT.BP(1))*(180/pi);
        
        % Find optimal segment orientation
        % -----------------------------------------------------------------
        PT.E_q = []; PT.E_e = []; PT.E_eA = []; PT.E_eB = [];
        options = optimset(fmin_settings{:}); % fminsearch options
        [~,exitval,exitflag] = fminsearch(@LOCAL_FitWorm,q0,options);
        
        % Append new fit data
        % -----------------------------------------------------------------
        PT.wx = [PT.wx PT.wx(end)+PT.x(end)];
        PT.wy = [PT.wy PT.wy(end)+PT.y(end)];
        PT.wz = [PT.wz PT.wz(end)+PT.z(end)];
        
        % Check if done
        % -----------------------------------------------------------------
        % Stop if all points fitted
        if sum(CHK_A) == length(CHK_A) && sum(CHK_B) == length(CHK_B),
            flag = false;
        end
        
        % Stop if error from single projection is too large
        eA = sqrt( (PT.Axc(end)-AAX).^2 + (PT.Ayc(end)-AAY).^2 );
        eB = sqrt( (PT.Bxc(end)-BBX).^2 + (PT.Byc(end)-BBY).^2 );
        if min(eA) > ER_thresh,
            disp(['eA error (',num2str(min(eA)),') is too large2'])
            flag = false;
            clip = find(PT.E_eA(:,end) <= ER_thresh,1,'last');
            %    JAE addition 140925
            wx = PT.wx;
            wy = PT.wy;
            wz = PT.wz;
            if isempty(clip), clip = -999; end
        elseif min(eB) > ER_thresh,
            disp(['eB error (',num2str(min(eB)),') is too large2'])
            flag = false;
            clip = find(PT.E_eB(:,end) <= ER_thresh,1,'last');
            %    JAE addition 140925
            wx = PT.wx;
            wy = PT.wy;
            wz = PT.wz;
            if isempty(clip), clip = -999; end
        end
        
        % Stop if fit folds back on itself
        dA = sqrt( (PT.Axc(end)-PT.Axc(1:end-PT.N)).^2 + ...
            (PT.Ayc(end)-PT.Ayc(1:end-PT.N)).^2 );
        dB = sqrt( (PT.Bxc(end)-PT.Bxc(1:end-PT.N)).^2 + ...
            (PT.Byc(end)-PT.Byc(1:end-PT.N)).^2 );
        % Commented out by @JAE 2014_03_06
        if min(dA) < ER_thresh && min(dB) < ER_thresh,
            disp('fit folding back on itself')
            flag = false;
            clip = 1;
            %    JAE addition 140925
            wx = PT.wx;
            wy = PT.wy;
            wz = PT.wz;
        end
        % Stop if output of fminsearch is NaN
        if isnan(exitval)
            fprintf('fminsearch returned a NaN, exiting fit\n')
            flag = false;
            wx = PT.wx;
            wy = PT.wy;
            wz = PT.wz;
        end
        
        % Clip bad data
        if isempty(clip), clip = 1; end
        if clip == -999,
            % Enters if no additional data is good to add -> clip it all
            if PT.N >= 2,
                PT.wx = PT.wx(1:(end-1)); % 1/3/2012
                PT.wy = PT.wy(1:(end-1)); % 1/3/2012
                PT.wz = PT.wz(1:(end-1)); % 1/3/2012
                % ---
                PT.Axc = PT.Axc(1:(end-(PT.N)));
                PT.Ayc = PT.Ayc(1:(end-(PT.N)));
                PT.Bxc = PT.Bxc(1:(end-(PT.N)));
                PT.Byc = PT.Byc(1:(end-(PT.N)));
            else
                PT.wx = PT.wx(1:(end-1));
                PT.wy = PT.wy(1:(end-1));
                PT.wz = PT.wz(1:(end-1));
                % ---
                PT.Axc = PT.Axc(1:(end-1));
                PT.Ayc = PT.Ayc(1:(end-1));
                PT.Bxc = PT.Bxc(1:(end-1));
                PT.Byc = PT.Byc(1:(end-1));
            end
        elseif ~isnan(clip) && PT.N >= 2,
            PT.wx(end) = PT.wx(end-1)+PT.x(clip);
            PT.wy(end) = PT.wy(end-1)+PT.y(clip);
            PT.wz(end) = PT.wz(end-1)+PT.z(clip);
            % ---
            PT.Axc = PT.Axc(1:(end-(PT.N-clip)));
            PT.Ayc = PT.Ayc(1:(end-(PT.N-clip)));
            PT.Bxc = PT.Bxc(1:(end-(PT.N-clip)));
            PT.Byc = PT.Byc(1:(end-(PT.N-clip)));
        elseif ~isnan(clip) && PT.N == 1,
            PT.wx = PT.wx(1:end-1);
            PT.wy = PT.wy(1:end-1);
            PT.wz = PT.wz(1:end-1);
            % ---
            PT.Axc = PT.Axc(1:(end-1));
            PT.Ayc = PT.Ayc(1:(end-1));
            PT.Bxc = PT.Bxc(1:(end-1));
            PT.Byc = PT.Byc(1:(end-1));
        end
        
        switch flip
            case 1, flip = 2;
            case 2, flip = 0;
        end
        
    end % if ~TGL_previousfit
    
    % Plot current step
    % ---------------------------------------------------------------------
    if TGL_PltSteps,
        figure(101); clf(101);
        set(101,'Position',[20 250 800 400])
        % ---
        subplot(1,2,1);
        plot(AAX,AAY,'k.'); hold on;
        plot(PT.Axc,PT.Ayc,'r*-','LineWidth',2);
        axis equal;
        title('FRONT View');
        % ---
        subplot(1,2,2);
        plot(BBX,BBY,'k.'); hold on;
        plot(PT.Bxc,PT.Byc,'r*-','LineWidth',2);
        axis equal;
        title('TOP View');
    end
    
end % while flag


% JAE addition: save the PT file so the two-dim back projections from the
% 3D-guessed points can be used for splining (useful if there is a
% basepoint mismatch)
PT.exitflag = exitflag;
PT.exitval = exitval;
%save 'PT_file' PT

%% Final plot
if TGL_PltFinal,
    figure;
    set(gcf,'Position',[20 250 800 400]);
    set(gcf,'Name','3D Merged Object','NumberTitle','off');
    % ---
    subplot(1,2,1);
    plot(AAX,AAY,'k.'); hold on;
    plot(PT.Axc,PT.Ayc,'r*-','LineWidth',2);
    legend({'Projection';'3D fit'}, ...
        'Location','SouthWest','FontSize',8);
    axis equal;
    title('FRONT View');
    % ---
    subplot(1,2,2);
    plot(BBX,BBY,'k.'); hold on;
    plot(PT.Bxc,PT.Byc,'r*-','LineWidth',2);
    axis equal;
    title('TOP View');
end

%% Return
wx = PT.wx;
wy = PT.wy;
wz = PT.wz;
PT.exitval = exitval;
PT.exitflag = exitflag;

%% Smoothed output
% -> Not yet tested
if TGL_spline,
    % Construct spline
    W = [wx;wy;wz];
    sp = spaps(1:length(wx),W,0);
    
    % Resample
    t = 1:1:length(wx);
    Wt = ppval(sp,t);
    
    % Plot
    figure;
    plot3(Wt(1,:),Wt(2,:),Wt(3,:),'b.-'); hold on;
    plot3(W(1,:),W(2,:),W(3,:),'ro');
    
    % Output
    wx = Wt(1,:)';
    wy = Wt(2,:)';
    wz = Wt(3,:)';
    
end


function e = LOCAL_FitWorm(q)
%% function e = LOCAL_FitWorm(q)
global DS PT mAAX mAAY mBBX mBBY A_proj B_proj A2B

% Target data
% -------------------------------------------------------------------------
% Find 3D point
dDS = (1/PT.N:1/PT.N:1).*DS;
dDY = zeros(1,length(dDS));
[x,y,z] = Get_RotateTranslate(dDS,dDY,[],[0 q(1) q(2)],[]);
PT.x = x; PT.y = y; PT.z = z;
x = [PT.wx PT.wx(end)+x];
y = [PT.wy PT.wy(end)+y];
z = [PT.wz PT.wz(end)+z];

% Compute reference projections: A_proj
[Axc,Ayc] = Get_3DtoCameraProjection(x,y,z,'proj',A_proj);

% Convert (x,y,z) from A_frame to B_frame
% ->  Y = rigid_motion(X,om,T)
r = rigid_motion([x;y;z],A2B{1},A2B{2});

% Compute reference projections; B_proj
[Bxc,Byc] = Get_3DtoCameraProjection(r(1,:),r(2,:),r(3,:),'proj',B_proj);

% Reference output
PT.Axc = Axc; PT.Ayc = Ayc;
PT.BP_A = [PT.Axc(end); PT.Ayc(end)];
PT.Bxc = Bxc; PT.Byc = Byc;
PT.BP_B = [PT.Bxc(end); PT.Byc(end)];

% Reference data
% -------------------------------------------------------------------------

% Plot check
if 0,
    figure(1000); clf(1000);
    subplot(1,2,1);
    plot(mAAX(1,:),mAAY(1,:),'k.'); hold on;
    plot(Axc(end-PT.N),Ayc(end-PT.N),'r*');
    plot(Axc,Ayc,'c.-');
    %    [~,temp_x,temp_y] = Draw_Circle(Axc(end-PT.N),Ayc(end-PT.N),PT.RR,'r');
    %    set(gca,'XLim',[min(temp_x) max(temp_x)]);
    %    set(gca,'YLim',[min(temp_y) max(temp_y)]);
    axis equal;
    
    subplot(1,2,2);
    plot(mBBX(1,:),mBBY(1,:),'k.'); hold on;
    plot(Bxc(end-PT.N),Byc(end-PT.N),'r*');
    plot(Bxc,Byc,'c.-');
    %    [~,temp_x,temp_y] = Draw_Circle(Bxc(end-PT.N),Byc(end-PT.N),PT.RR,'r');
    %    set(gca,'XLim',[min(temp_x) max(temp_x)]);
    %    set(gca,'YLim',[min(temp_y) max(temp_y)]);
    axis equal;
end

% Compute errors
% -------------------------------------------------------------------------

% Compute error for view A
E = sqrt( ...
    (repmat(Axc(end-PT.N+1:end)',1,size(mAAX,2))-mAAX).^2 + ...
    (repmat(Ayc(end-PT.N+1:end)',1,size(mAAX,2))-mAAY).^2);
eA = min(E,[],2);

% Compute error for view B
E = sqrt( ...
    (repmat(Bxc(end-PT.N+1:end)',1,size(mBBX,2))-mBBX).^2 + ...
    (repmat(Byc(end-PT.N+1:end)',1,size(mBBX,2))-mBBY).^2);
eB = min(E,[],2);

% Sum errors to get final error
e = sum(sum(eA) + sum(eB));

% Error metric for analysis
PT.E_q = [PT.E_q; q];
PT.E_e = [PT.E_e; e];
PT.E_eA = [PT.E_eA, eA];
PT.E_eB = [PT.E_eB, eB];


%% function V = LOCAL_ClockGuess(Rx,Ry,x,y,RR)
function V = LOCAL_ClockGuess(Rx,Ry,x,y,RR)
thresh = RR/5;
alpha = 30; % degrees of wedge to cut data

% Correct orientation
if size(Rx,1) == 1, Rx = Rx'; end
if size(Ry,1) == 1, Ry = Ry'; end

% Upsample and spread out (x,y)
t = 0:1/(length(x)-1):1;
tt = -0.5:0.01:1;
xx = interp1(t,x,tt,'linear','extrap');
yy = interp1(t,y,tt,'linear','extrap');

[xxt,yyt] = rotate2(xx,yy,+alpha*(pi/180),[x(end) y(end)]);
[xxb,yyb] = rotate2(xx,yy,-alpha*(pi/180),[x(end) y(end)]);
[xxtt,yytt] = rotate2(xx,yy,+alpha/2*(pi/180),[x(end) y(end)]);
[xxbb,yybb] = rotate2(xx,yy,-alpha/2*(pi/180),[x(end) y(end)]);

xx = [xx,xxt,xxtt,xxb,xxbb];
yy = [yy,yyt,yytt,yyb,yybb];

% Find distances
XX = repmat(xx,length(Rx),1);
YY = repmat(yy,length(Ry),1);
RRx = repmat(Rx,1,length(xx));
RRy = repmat(Ry,1,length(xx));
D = sqrt( (XX-RRx).^2 + (YY-RRy).^2 );
D = min(D,[],2);

% Threshold
good = logical( D > thresh);
Rxx = Rx(good);
Ryy = Ry(good);

% Check if data exists
if isempty(Rxx),
    V = [x(end);y(end)];
    return;
end

% Fit slope to data
m = polyfit(Rxx,Ryy,1);
th_guess = atan2(m(1),1);

% Find optimal theta
C = [x(end),y(end)];
options = optimset('maxfunevals',500,'MaxIter',500); % fminsearch options
[th,exitval,exitflag] = fminsearch(@LOCAL_FitClock,th_guess,options,Rxx,Ryy,RR,C);

% Output final V
V = [0;0];
V(1) = RR*cos(th) + C(1);
V(2) = RR*sin(th) + C(2);

% Plot check
if 0,
    figure(1001); clf(1001);
    plot(Rx,Ry,'k.'); hold on;
    plot(Rxx,Ryy,'b.');
    plot(xx,yy,'c.');
    plot(x,y,'m*-');
    plot([C(1) V(1)],[C(2) V(2)],'r*-','LineWidth',2);
    Draw_Circle(C(1),C(2),RR,'r');
    axis equal;
end
%% function e = LOCAL_FitClock(q,Rx,Ry,R,C)
function e = LOCAL_FitClock(q,Rx,Ry,R,C)
x(1) = C(1);
y(1) = C(2);
x(2) = R*cos(q) + C(1);
y(2) = R*sin(q) + C(2);

% Threshold distance (Rx,Ry) to line
% -> see http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
% -> Eqn #14
%       (x2,y2) = (x(2),y(2)
%       (x1,y1) = (x(1),y(1))
%       (x0,y0) = (Rx,Ry)
% d = abs( (x(2)-x(1)).*(y(1)-Ry) - (x(1)-Rx).*(y(2)-y(1)) )./ ...
%     sqrt( (x(2)-x(1))^2 + (y(2)-y(1))^2);
% e = sum(d);

d = abs( (Rx-x(2)).^2 + (Ry-y(2)).^2 );
e = sum(d);

if 0,
    figure(1002); clf(1002);
    plot(Rx,Ry,'k.'); hold on;
    plot(C(1),C(2),'c*');
    plot(x,y,'r*-');
    Draw_Circle(C(1),C(2),R,'r');
end