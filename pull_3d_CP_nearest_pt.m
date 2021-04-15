function CP = pull_3d_CP_nearest_pt(PTfile,threeDwhisker,top_down_contact,contact_binary);
%% function CP = pull_3d_CP_nearest_pt(PTfile,3Dwhisker,top_down_contact)
%
%   Method:
%       Find the two closest points on the three-dimensional whisker to the
%       contact point in the top-down image and return the point on the
%       line between those two which corresponds to the x-coordinate of the
%       top-down contact point
%
%   Input:
%       top_down_contact    = n x 2 array of contact points from the top-down
%                               camera
%       PTfile              = the 'summary_PT' file saved from Merge3D
%                               (Ellis edits provide this output, otherwise save it from within the command script)
%       threeDwhisker       = tracked whisker in three dimensions 
%       contact_binary      = logical vector indicating contact periods
%

for ii = 1:length(PTfile)
    
    if contact_binary(ii)
        for pt = 1:length(PTfile{ii}.Bxc)

            distance_to_topCP(pt) = abs(PTfile{ii}.Bxc(pt) - top_down_contact(ii,1));
        end
        
        % Grab the segmentof 3D whisker closest to our point
        closest_x = PTfile{ii}.Bxc(distance_to_topCP==min(distance_to_topCP));
        closest_y = PTfile{ii}.Byc(distance_to_topCP==min(distance_to_topCP));
        distance_to_topCP(distance_to_topCP==min(distance_to_topCP)) = [];
        next_closest_x = PTfile{ii}.Bxc(distance_to_topCP==min(distance_to_topCP));
        next_closest_y = PTfile{ii}.Byc(distance_to_topCP==min(distance_to_topCP));
        
        pfit = polyfit([closest_x,next_closest_x],[closest_y,next_closest_y],1);
        CP{ii} = polyval(pfit,top_down_contact(ii,1));
        
%         if min(distance_to_topCP) == min(distance_to_topCP(1:length(threeDwhisker{ii})))
%             CP{ii} = threeDwhisker{ii}(:,distance_to_topCP==min(distance_to_topCP(1:length(threeDwhisker{ii}))));
%         else
%             bads(ii) = 1;
%         end
        
    end
    
end

disp(['Frames where the contact point is not in tracked 3D whisker: ',num2str(sum(bads))])