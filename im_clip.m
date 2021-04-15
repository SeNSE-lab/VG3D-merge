function q = im_clip(q0,v,dir,varargin)
%% function q = im_clip(q0,v,dir,{'param_name',param})
% -------------------------------------------------------------------------
% INPUT:
%   q0 - original image
%   v - (x,y) coordinate of pivot point, format: [x;y]
%   dir - quadrant to clip, string of numbers '1','2','3','4','12',etc ...
%   varargin:
%       'plot' - toggle plot
%       'fv' - fill value, {default: 0}
% OUTPUT:
%   q - clipped image
% -------------------------------------------------------------------------
% Brian Quist
% November 10, 2011

%% Handle inputs
fv = 0;
TGL_plot = 0;

% User inputs
% -------------------------------------------------------------------------
if ~isempty(varargin),
    for ii = 1:2:(length(varargin))
        switch varargin{ii},
            case 'fv',      fv = varargin{ii+1};
            case 'plot',    TGL_plot = varargin{ii+1};
            otherwise,
                error('Not a valid input parameter');
        end
    end
end

%% Cycle through quadrants
q = q0;
for ii = 1:length(dir)
    switch dir(ii)
        case '1', q(1:v(2),v(1):end) = fv; 
        case '2', q(1:v(2),1:v(1)) = fv;
        case '3', q(v(2):end,1:v(1)) = fv; 
        case '4', q(v(2):end,v(1):end) = fv;
        otherwise
            error('Not a valid quadrant');
    end
end

%% Check plot
if TGL_plot,
    figure;    
    imshow(q); hold on;
    plot(v(1),v(2),'r*');
    plot([1 size(q,2)],[v(2) v(2)],'r');
    plot([v(1) v(1)],[1 size(q,1)],'r');
end