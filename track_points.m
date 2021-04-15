function trackedpoints=track_points(frame,step_btwn_pts, perp_line_lgth, startx,starty)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Written by Matthew Graff
%   Last Modified 1/31/12
%
%   This function takes an image with a single whisker in it and tracks 
%    its position. For best results, the images should have a black whisker
%     on a white background. The higher the contrast the better.
%
%
%   INPUTS:
%   frame- an image matrix of a whisker running left to right
%
%
%   step_btwn_pts- This is the approximate number of pixels that you want
%       between points on the whiskers.
%
%   prep_line_lgth- This is the length of the perpendicular line which
%       samples the intensities across the whisker to find the darkest
%       pixel.
%
%   startx- x-coord of where you want the tracking code to start
%
%   starty- y-coord of where you want the tracking code to start
%
%   OUTPUTS:
%    trackedpoints- an array of x-y tracked values along the whisker
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% warning off;

% change plot_tracking to 1 for debugging
% frame=histeq(frame);


plot_tracking=1;

num_pts=round((size(frame,1))/step_btwn_pts);

trackedpoints=zeros(num_pts,2);
pixvals=zeros(perp_line_lgth*2+1,num_pts);


x=round(startx);
y=round(starty);
slope=0;

for i=1:num_pts
    
    
    
%     if plot_tracking
%         if i==1
%             figure;
%             imshow(frame);
%         else
%             hold on;
%             impoint(gca, x, y);
%             c=imline(gca, [x,x], [y-perp_line_lgth,y+perp_line_lgth]);
%             
%             %     pause; %uncomment to see point by point tracking or add breakpoint
%         end
%     end
    
    % newer quick code
    if (x<=1024&&x>=1)&&(y<=1024&&y>=1)
        if slope==0;
            slope2=-999999;
        end
        if slope2<0
            sign=-1;
        else
            sign=1;
        end
        LL=2*perp_line_lgth+1;
        dx=1/sqrt(1+slope2^2);
        dy=sign*sqrt(1^2-dx^2);
        mx=x-perp_line_lgth*dx:dx:x+perp_line_lgth*dx;
        my=y-perp_line_lgth*dy:dy:y+perp_line_lgth*dy;
        valsx=round(mx);
        valsy=round(my);
        for j=1:LL
            if valsx(j)>1024
                valsx(j)=1024;
            elseif valsx(j)<1
                valsx(j)=1;
            end
            if valsy(j)>1024
                valsy(j)=1024;
            elseif valsy(j)<1
                valsy(j)=1;
            end
            pixvals(j,i)=frame(valsy(j),valsx(j));
        end
        

        
        [~,ind]=min(pixvals(:,i));
        
        
        trackedpoints(i,1)=valsx(ind);
        trackedpoints(i,2)=valsy(ind);
        
        
        
        
        
        if i>10
            if trackedpoints(i,1)-trackedpoints((i-10),1)==0
                slope = 999999;
            else
                slope = (trackedpoints(i,2)-trackedpoints((i-10),2))...
                    /(trackedpoints(i,1)-trackedpoints((i-10),1));
            end
            
            slope2 = -1/slope;
        end
        
%         if plot_tracking
%             d=impoint(gca,x, valsy(ind));
%             setColor(d,'red');
%         end
        
        x=round(valsx(ind)+sqrt(step_btwn_pts^2/(slope.^2 +1)));
        y=round(valsy(ind)+slope*sqrt(step_btwn_pts^2/(slope.^2 +1)));
        
        
        
    else
        break;
    end
end
trackedpoints=trackedpoints(1:i-1,:);


end
