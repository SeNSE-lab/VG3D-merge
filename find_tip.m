function [tracks_trunc, int_vals]=find_tip(frame, trackedpoints,thresh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by: Matthew Graff
% Last modified: 1/31/12
%
%  Inputs-
%      frame-   The image from which the tracked points came from
%               trackedpoints- nx2 matrix of tracked points (ouput of
%               track_points.m function)
%
%      thresh-  the intensity threshold for which indicates that
%               the tracking code has left the whisker
%
%      est_tip- Estimated index for the tip of the whisker. This variable
%               should be the lower bound for where the code should start
%               looking for the tip. This is helpful if there are any
%               strange peaks in intensity before the tip. Also speeds up
%               the code slightly.
%
%      extend-  if find_tip outputs tracks that are consistantly
%               a few pixels too long/short this variable will add or
%               remove additional points. Pos num to add points Neg
%               num to remove.
%
% Outputs-
%     tracks_trunc- The truncated version of the input trackedpoints
%                   where the last point is the whiskers tip
%
%     int-vals-     the intensity values of the tracked pixels. Plot
%                   this to get an idea of how to set thresh.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


image_corrected=histeq(frame);
round(trackedpoints);
int_vals=zeros(1,size(trackedpoints,1));
for i=1:size(trackedpoints,1)
    int_vals(i)=image_corrected(trackedpoints(i,2),trackedpoints(i,1));
end
i=find(int_vals(1,1:end)>thresh,1,'first');
% i=est_tip+i+extend;
tracks_trunc=trackedpoints(1:i,:,:);
end