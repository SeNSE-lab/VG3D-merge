function [new_xpts,new_ypts] = check_whisker_smoothness(xpts,ypts,num_stdevs,threshold)
%%  function check_whisker_smoothness()
%   Checks to make sure there is no unbecoming weirdness in whiskers due
%   to poor tracking, and smooths if necessary

% Determine number of good points

diffy = diff(ypts);
good_points = logical(diffy < (mean(diffy)+std(diffy)*num_stdevs) & diffy > (mean(diffy)-std(diffy)*num_stdevs));

% Base Case: good whisker with few bad points
if length(xpts) - sum(good_points) <= threshold
    new_xpts = xpts;
    new_ypts = ypts;
    
else
    good_points(1)=1;good_points(end+1)=1;
    good_xpts = xpts(good_points);
    good_ypts = ypts(good_points);
    pfit = polyfit(good_xpts,good_ypts,2);
    new_xpts = xpts(1):xpts(end);
    new_ypts = polyval(pfit,new_xpts);
end