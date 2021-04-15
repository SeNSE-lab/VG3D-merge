function q1 = im_cleanborder(q0,varargin)
%% function q1 = im_cleanborder(q0,{'param_name',param})
% -------------------------------------------------------------------------
% INPUT:
%   q0 - input image
%   varargin:
%       'bwidth' - border width, in pixels, to remove {default: 2 pixels}
%       'fill_v' - value to fill border {default: 0}
%       'plot' - plot results {default: false}
% -------------------------------------------------------------------------
% Brian Quist
% October 20, 2011

%% Handle inputs
bwidth = 2;
fill_v = 0; 
TGL_plot = 0; 
if ~isempty(varargin),
   for ii = 1:2:(length(varargin)) 
       switch varargin{ii},
           case 'bwidth', bwidth = varargin{ii+1};
           case 'fill_v', fill_v = varargin{ii+1};
           case 'plot', TGL_plot = varargin{ii+1};
           otherwise,
               error('Not a valid input parameter');
       end
   end
end

%% Input check
bwidth = round(bwidth);
if bwidth < 1, bwidth = 1; end

%% Fill values:
q1 = q0;
[a,b] = size(q0);
q1(:,1:bwidth) = fill_v;
q1(:,(b-bwidth+1):b) = fill_v;
q1(1:bwidth,:) = fill_v;
q1((a-bwidth+1):a,:) = fill_v;

%% Plot
if TGL_plot,
   figure;
   subplot(1,2,1);
   imshow(q0);
   title('Original');
   subplot(1,2,2);
   imshow(q1);
   title('Processed');    
end
