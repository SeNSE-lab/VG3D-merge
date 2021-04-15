function CPf = get_front_contact_pt(PTfile,peg_horiz_ln,contact_frames);
%%  function CPf = get_front_contact_pt(whisker_frames,peg_horiz_ln)
%   uses a single x-coordinate to look for y-coords of contact using a
%   <logical> contact frame list, usually from the top-down view

for ii = 1:length(PTfile)
    xf = PTfile{ii}.Axc;
    yf = PTfile{ii}.Ayc;
    
    % If in contact, extract x_coord
    if contact_frames(ii)
        old_distance = Inf;
        for link = 1:length(xf)
            new_distance = abs(306 - xf(link));
            if old_distance < new_distance
                x_coord = xf(link-1);
                break
            else
                old_distance = new_distance;
            end
        end
%         min_x_diff = min(abs(x_coords));
%         for pp = 1: length(x_coords)
%             if min(abs(x_coords(pp))) == min_x_diff
%                 x_pt = x_coords(pp);
%                 break
%             end
%         end
       
        CPf(ii,:) = [x_coord,yf(xf==x_coord)];
    else
        CPf(ii,:) = [NaN,NaN];
    end
end
