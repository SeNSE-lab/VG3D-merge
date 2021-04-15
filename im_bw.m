function q1 = im_bw(q0,varargin)
%% function q1 = im_bw(q0,{'setting_name',setting})
% -------------------------------------------------------------------------
% INPUT:
%   q0 - original gray-scale image
%   varargin:
%       'invert' - invert bw selection after threshold
%           -> {Default: false}
%       'bw_thresh' - bw threshold
%           -> {Default: 0.2}
%       'max_prop' - only return maximum of certain regionprops type
%           -> {Default: 'none'};
%           -> any regionprops name is valid
%       'max_prop_n' - number of regionprop areas to include in return image
%           -> {Default: 1}
%       'plot' - toggle output plot
%           -> {Default: 0} or 1
% OUTPUT:
%   q1 - bw thresholded image
% -------------------------------------------------------------------------
% Brian Quist
% November 10, 2011

%% Handle inputs
bw_thresh = 0.2;
max_prop = 'none';
max_prop_n = 1;
TGL_plot = 0; 
TGL_invert = 0;
if ~isempty(varargin),
   for ii = 1:2:(length(varargin)) 
       switch varargin{ii},
           case 'bw_thresh', bw_thresh = varargin{ii+1};
           case 'max_prop', max_prop = varargin{ii+1};
           case 'max_prop_n', max_prop_n = varargin{ii+1};
           case 'plot', TGL_plot = varargin{ii+1};
           case 'invert', TGL_invert = varargin{ii+1};
           otherwise,
               error('Not a valid input parameter');
       end
   end
end

%% Threshold image
q1 = im2bw(q0,bw_thresh);

if TGL_invert,
    q1 = ~q1;
end

%% Select out regionprop
if ~strcmp('none',max_prop),
    
    % Run region_props
    lbl = bwlabel(q1);
    rp = regionprops(lbl,max_prop);
    
    % Sort
    rp = struct2array(rp);
    rp = [rp;1:length(rp)]';
    rp = sortrows(rp,-1);
    
    % Filter regions
    q1 = zeros(size(q1));
    for ii = 1:max_prop_n,
        if max_prop_n <= size(rp,1),
            keep = logical(lbl == rp(ii,2));
            q1(keep) = 1;
        end
    end
    
    % Convert back to binary
    q1 = logical(q1);
end

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
