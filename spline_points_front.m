function [xf,yf] = spline_points_front(xpts,ypts)
%% Must have the points (x and y) loaded into variable space to proceed

if iscell(xpts)
    for ii = 1:length(xpts)
        notNaNs = xpts{ii}(isnan(xpts{ii})==0);
        xf{ii} = [min(notNaNs):max(notNaNs)];
        yf{ii} = spline(xpts{ii},ypts{ii},xf{ii}');
    end
else
    for ii = 1:size(xpts,2)
        notNaNs = xpts(isnan(xpts(:,ii))==0,ii);
        xf{ii} = [min(notNaNs):max(notNaNs)];
        yf{ii} = spline(xpts(:,ii),ypts(:,ii),xf{ii}');
    end
end