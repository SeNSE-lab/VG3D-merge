function tracked_3D = mergeLoop(f,t,calibInfo,tracked_3D_fileName,varargin)
%% function mergeLoop(f,t,calibInfo,tracked_3D_fileName,[count],[onlyC],[plotTGL],[C])

%
% INPUTS:
%             t = top whisker struct
%             f = front whisker struct
%
%             calibInfo =1x10 cell of calibration information
%
%             [count] = if the data have already been partially merged,
%               tells the function where to start.
%             [onlyC] = Logical to indicate if we should merge only the
%               contact frames. Defaults to 1.
%             [plotTGL] = Logical that switches to a viewing only mode. If
%               plotTGL is true, no data is saved. Default is 0.
%             [C] = A contact binary if supplied. Defaults to all nans
%
%
% Saves tracked 3D whisker struct to a file in the current directory with tag
% 'tracked_3D'. e.g. :'rat2015_15_JUN11_VG_D4_t01_tracked_3D.mat'


%% Input handling
numvargs = length(varargin);
narginchk(4,8)
optargs = {1,true,false,nan(length(f),1)};

% overwrite user supplied args
optargs(1:numvargs) = varargin;
[count,onlyC,plotTGL,C] = optargs{:};

%% Output file naming
if exist(tracked_3D_fileName,'file')
    load(tracked_3D_fileName);
end


%% Set defaults and input errors

assert(length(f) == length(t));
assert(length(f) == length(C));

notify = 0;
count_at_error = 0;
numErrors = 0;
numFrames = length(f);
if ~exist('tracked_3D','var')
    tracked_3D = struct([]);
end
step = 1000;% size of par loop
save_step = 1; % save after every 11 clips
DS = 1;% sets the minimum internode distance.
minWhiskerSize = 20; % in # of nodes
N = 20; % More should give a stabler fit.
ext = 5; % how many frames to merge before and after contact
merging = true;
clips = 1:step:numFrames;

if onlyC
    mergeFlags = LOCAL_getMergeFlags(C,numFrames,ext);
else
    % merge every frame
    mergeFlags = true(size(C));
end

%%
disp('merging...')
% Outer loop is big serial chunks that saves every [step] frames
while merging
    currentClip = clips(count);
    % Makes sure we don't try to access a frame past the last frame.
    if currentClip == clips(end)
        iter = numFrames-clips(end);
    else
        iter = step-1;
    end
    % if plotting is requested, use this loop. It does not save anything,
    % nor actually compute. The primary goal is fo make sure the
    % calibration is correct.
    if plotTGL
        for i = currentClip:currentClip+iter
            if ~mergeFlags(i) || isempty(t(i).x) || isempty(f(i).x)
                tracked_3D(i).x = []; tracked_3D(i).y = []; tracked_3D(i).z = [];
                tracked_3D(i).time = i-1;
                continue
            end
            
            [merge_x,merge_y,merge_z]= Merge3D_JAEv1(f(i).x,f(i).y,t(i).x,t(i).y,i,calibInfo,'plot',1,'wm_opts',{'DS',DS,'N',N});
        end
        close all
        continue
    end
    
    % Parallel for loop which does the actual merging. Gets batches from
    % the current outer loop.
    try
        parfor i = currentClip:currentClip+iter
            fprintf('Merging frame \t%i of \t%i\n',i,numFrames)
            %initialize the merged values in the parfor loop.
            merge_x = [];merge_y = [];merge_z = [];last_merge_x = []; last_merge_y = []; last_merge_z = [];
            
            % if this frame is not flagged for merging, skip it. This should
            % have already checked for empties in both views
            
            if ~mergeFlags(i) || isempty(t(i).x) || isempty(f(i).x)
                tracked_3D(i).x = []; tracked_3D(i).y = []; tracked_3D(i).z = [];
                tracked_3D(i).time = i-1;
                continue
            end
            
            
            % Initial merge.
            [merge_x,merge_y,merge_z]= Merge3D_JAEv1(f(i).x,f(i).y,t(i).x,t(i).y,i,calibInfo,'plot',0,'wm_opts',{'DS',DS,'N',N});
            
            x_out = merge_x;
            y_out = merge_y;
            z_out = merge_z;
            
            tracked_3D(i).x = x_out; tracked_3D(i).y = y_out; tracked_3D(i).z = z_out;
            tracked_3D(i).time = i-1;tracked_3D(i).frontTime = f(i).time;tracked_3D(i).topTime = t(i).time;
            
            close all
        end
        count = count + 1;
    catch
        notify = 1;
        count_at_error = count;
        numErrors = numErrors+1;
        if numErrors == 1;
	    disp('Error, attempting restart')
        elseif numErrors > 1
	    disp('There have been repeated errors in merging. Skipping segment')
            for i = currentClip:currentClip+iter
                tracked_3D(i).x = []; tracked_3D(i).y = []; tracked_3D(i).z = [];
                tracked_3D(i).time = i-1;tracked_3D(i).frontTime = [];tracked_3D(i).topTime = [];
            end
            count = count +1;
            notify = 0;
            numErrors = 0;
            
        end
        gcp
    end
    if count > count_at_error & notify
        disp('GCP restart successful')
        notify = 0;
        numErrors = 0;
    end
    
    if count>length(clips)
        merging = 0;
    end
    if mod(count,save_step) == 0
        save(tracked_3D_fileName,'tracked_3D','C','calibInfo','count')
    end
end
save(tracked_3D_fileName,'tracked_3D','C','calibInfo','count')
end



%% Local function to pad the contact variable by a certain number of frames.
% This allows us to only merge a little before and after contact.

function mergeFlags = LOCAL_getMergeFlags(C,numFrames,ext)
%% function mergeFlags = LOCAL_getMergeFlags(C,numFrames,ext)


mergeFlags = C;
cc = convertContact(C);

for ii = 1:size(cc,1)
    if (cc(ii,1)-ext)<1
        mergeFlags(1:cc(ii,2)+ext) = 1;
    elseif (cc(ii,2)+ext)>numFrames
        mergeFlags(cc(ii,1)+ext:end);
    else
        mergeFlags(cc(ii,1)+ext:cc(ii,2)+ext);
    end
end

end

