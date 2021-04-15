function [new_xpts,new_ypts] = extend_top_down_whisker(xpts,ypts,whisker_side,extent)
%%  function [new_xpts,new_ypts] = extend_top_down_whisker(xpts,ypts,extend_dir,extent)
%   Takes a whisker from the top-down camera view and extends it towards
%   the base
%
%   Input:  xpts, ypts  = vectors
%               extent  = pixel number to extend
%
% 	IMPORTANT!  Whisker should be input with the xpts(1) < xpts(end);
% 	increasing from first to last entry on the x-axis

pfit = polyfit(xpts,ypts,2);

switch whisker_side
    case 'L'
        new_xpts = [xpts(1)-extent:xpts(end)];
        new_ypts = polyval(pfit,new_xpts);
        
    case 'R'
        new_xpts = [xpts(1):xpts(end)+extent];
        new_ypts = polyval(pfit,new_xpts);
        
end