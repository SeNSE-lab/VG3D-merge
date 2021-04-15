function [xc,yc] = Get_centerline(x,y,varargin)
%% function [xc,yc] = Get_centerline(x,y,{'param_name',param})
% -------------------------------------------------------------------------
% INPUT:
%   (x,y) - coordinates of line-object
%   varargin:
%       'N' - number of points to use in approximation
%       'R' - radius of circle used in estimating slope
%       'W' - width of segment to grab centerline approximation
%       'plot' - toggle to plot
% OUTPUT:
%   (xc,yc) - coordinates of centerline for line-object
% -------------------------------------------------------------------------
% NOTES:
%   + Requires: 
%       -> rotate2
% -------------------------------------------------------------------------
% Brian Quist
% November 1, 2011

%% Initial settings
TGL_plot = 0;
% --------------
N = round(length(x)*0.05);
R = sqrt((x(end)-x(1))^2+(y(end)-y(1))^2)*0.01;
W = R/10;
if N < 15 && length(x) >= 15, N = 15; end 

%% Handle inputs
if ~isempty(varargin),
   for ii = 1:2:(length(varargin)) 
       switch varargin{ii},
           case 'N', N = varargin{ii+1};
           case 'R', R = varargin{ii+1};
           case 'W', W = varargin{ii+1};
           case 'plot', TGL_plot = varargin{ii+1};
           otherwise,
               error('Not a valid input parameter');
       end
   end
end

%% Get approximate centerline
index = round(1:(length(x)-1)/(N-1):length(x));
xc0 = x(index);
yc0 = y(index);
xc = xc0;
yc = yc0;

%% Cycle through and adjust centerline
for ii = 1:length(xc),
    
    % Setup filter
    good = logical( sqrt((x-xc(ii)).^2 + (y-yc(ii)).^2) <= R );
    
    % Linear-fit data for trend
    p = polyfit(x(good),y(good),1);
    
    % Rotate local data and center and (xc,yc)
    th = atan2(p(1),1);
    xr = x(good)-xc(ii);
    yr = y(good)-yc(ii);
    [xr,yr] = rotate2(xr,yr,-th,[0 0]);    
    
    % Find max and min near region
    goodw = logical( -W <= xr & xr <= W );
    pt = mean( yr(goodw) );
    
    % Reposition
    [xn,yn] = rotate2(0,pt,th,[0 0]);
    xn = xn + xc(ii);
    yn = yn + yc(ii);
    
    % Replace
    xc(ii) = xn;
    yc(ii) = yn;
    
    % Plot check
    if 0,
       figure(101); clf(101);
       subplot(1,2,1);
       plot(x,y,'k.'); hold on;
       plot(x(good),y(good),'b.');
       plot(xc0,yc0,'r.-');
       plot(xc0(ii),yc0(ii),'r*');
       plot(xc(ii),yc(ii),'c*');
       subplot(1,2,2);
       plot(xr,yr,'k.'); hold on;
       plot(xr(goodw),yr(goodw),'r.');
       plot(0,pt,'c*'); 
       axis equal;
    end
    
end

%% Final plot check
if TGL_plot,
   figure;
   plot(x,y,'k.'); hold on;
   plot(xc,yc,'r.-');
end