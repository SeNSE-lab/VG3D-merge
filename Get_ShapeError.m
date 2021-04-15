function [e,stats] = Get_ShapeError(xa,ya,za,xb,yb,zb,varargin)
%% [e,stats] = Get_ShapeError(xa,ya,za,xb,yb,zb,{'param_name',param})
% -------------------------------------------------------------------------
% INPUT:
%   (xa,ya,za) - 3D coordinates of object A
%   (xb,yb,zb) - 3D coordinates of object B
%   * Any one input set to [] will be filled with zeros
%   varargin:
%       'plot' - toggle plotting
%       'tgl_e' - toggle which object to compute error w.r.t {default: 'B'}
%               --> either 'A' or 'B'
% OUTPUT:
%   e - error w.r.t. object {default: w.r.t. object B}
%   stats - struct of error fit stats
%       .max - maximum error
%       .min - minimum error
%       .median - median error
%       .std - standard deviation of error
%       .N - number of points in calculation
% -------------------------------------------------------------------------
% NOTES:
%   + Helper functions:
%       -> arclength3d
% -------------------------------------------------------------------------
% Brian Quist
% November 7, 2011

%% Handle variable inputs
TGL_plot = 1;
TGL_objE = 'B';
if ~isempty(varargin),
   for ii = 1:2:(length(varargin)) 
       switch varargin{ii},
           case 'tgl_e',    TGL_objE = varargin{ii+1};
           case 'plot',     TGL_plot = varargin{ii+1};
           otherwise,
               error('Not a valid input parameter');
       end
   end
end

%% Structure inputs
Na = max([length(xa),length(ya),length(za)]);
Nb = max([length(xb),length(yb),length(zb)]);

% Handle missing data
if isempty(xa), xa = zeros(1,Na); end
if isempty(ya), ya = zeros(1,Na); end
if isempty(za), za = zeros(1,Na); end
if isempty(xb), xb = zeros(1,Nb); end
if isempty(yb), yb = zeros(1,Nb); end
if isempty(zb), zb = zeros(1,Nb); end

% Correct vector size
% object A --> row vector
% object B --> column vector
if size(xa,2) == 1, xa = xa'; end
if size(ya,2) == 1, ya = ya'; end
if size(za,2) == 1, za = za'; end
if size(xb,1) == 1, xb = xb'; end
if size(yb,1) == 1, yb = yb'; end
if size(zb,1) == 1, zb = zb'; end


%% Create error calculation maxtix
XA = repmat(xa,Nb,1);
YA = repmat(ya,Nb,1);
ZA = repmat(za,Nb,1);
XB = repmat(xb,1,Na);
YB = repmat(yb,1,Na);
ZB = repmat(zb,1,Na);

%% Compute distances
D = sqrt( (XA-XB).^2 + (YA-YB).^2 + (ZA-ZB).^2 );

%% Compile error metrics and arclength
switch TGL_objE
    case 'A',
        dim = 1;
        [s_total,s] = arclength3d(xa,ya,za);
    case 'B',
        dim = 2;
        [s_total,s] = arclength3d(xb,yb,zb);
    otherwise,
        error('Bad tgl_e selection');
end

%% Output
e = min(D,[],dim);
stats.max = max(e);
stats.min = min(e);
stats.mean = mean(e);
stats.median = median(e);
stats.std = std(e);
stats.N = length(e);
stats.s = s;
stats.s_total = s_total;

%% Plot
if TGL_plot,
   figure;
   set(gcf,'Name','3D Error','NumberTitle','off');
   plot(s,e,'b.-'); hold on;
   plot([s(1) s(end)],[stats.max stats.max],'k:');
   plot([s(1) s(end)],[stats.mean stats.mean],'r-');
   plot([s(1) s(end)],[stats.mean+stats.std stats.mean+stats.std],'m:');
   plot([s(1) s(end)],[stats.mean-stats.std stats.mean-stats.std],'m:');
   xlabel(['Arc length along object ',TGL_objE]);
   ylabel('Error (abs)');
   text(0.01,0.99, ...
       {['Max(E): ',num2str(stats.max),' ... k']; ...
        ['Mean(E): ',num2str(stats.mean),' ... r']; ...
        ['Std(E): ',num2str(stats.std),' ... m']; ...
        ['--------']; ...
        ['Percent (of 100) of s):']; ...
        ['Max(E): ',num2str(stats.max/stats.s_total*100)]; ...
        ['Mean(E): ',num2str(stats.mean/stats.s_total*100)]; ...
        ['Std(E): ',num2str(stats.std/stats.s_total*100)]},...
        'Units','normalized','FontSize',8, ...
        'VerticalAlignment','top','BackgroundColor','w'); %#ok<NBRAK>
end