function [xpts,ypts] = bridge_gaps(xpts,ypts)
%% function bridge_gaps(xpts,ypts)
% bridge gaps in (usually the top view's) tracked whisker

for ii = 1:length(xpts)
    xpts{ii}=double(xpts{ii});
    for p = 2:length(xpts{ii})
        if abs(xpts{ii}(p)-xpts{ii}(p-1)) > 3 && ~(xpts{ii}(p)-xpts{ii}(p-1) == 0)
            if mean(diff(xpts{ii})) < 0.8
                ypts{ii}=ypts{ii}(1:2:end);
                xpts{ii}=xpts{ii}(1:2:end);
            end
            ypts{ii} = spline(xpts{ii},ypts{ii},xpts{ii}(1):xpts{ii}(end));
            xpts{ii} = xpts{ii}(1):xpts{ii}(end);
            break
        end
    end
end