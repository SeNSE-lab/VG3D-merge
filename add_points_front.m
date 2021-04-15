function [pfit_x,pfit_y] = add_points_front(xf,yf,image_span)
%%  Fix manually tracked whiskers that did not extend to tip
%   Have the xf and yf (e.g. x-front, y-front) data loaded 
%   Be IN the directory of the images
%   image_span should be a vector of image numbers to process
global ImgNum
global prefix
ImgNum = 1;

% Set the degree of the polyfit polynomial
polyfit_degree = 3;

filename = ls('*000001.tif');
prefix = filename(1:25);

fig = figure(10201); % call palandromic figure

for img = image_span
    thisnum = sprintf('%06d',img); % pad image number
    imshow(imread([prefix,thisnum,'.tif'])); % concatenate prefix and image number... print picture: rat
    title(sprintf('click on exta points\nPress ENTER to end\nPress BACKSPACE to delete previous point')) % tell the people what they need to know
    uicontrol('Style','text','Position',[400 50 60 20],'String',['Frame: ',num2str(img)]); % show image number
    hold on % stay!
    plot(xf(:,img),yf(:,img)) % plot the already-tracked part of the whisker
    [x_to_tip{img},y_to_tip{img}] = getpts(fig);
    
    % Polyfit the results
    splined_y{img} = spline([xf(:,img);x_to_tip{img}],[yf(:,img);y_to_tip{img}],min(xf(isnan(xf(:,img))==0,img)):x_to_tip{img}(end));
    pfit_x{img} = min(xf(isnan(xf(:,img))==0,img)):x_to_tip{img}(end);
    pfit = polyfit(min(xf(isnan(xf(:,img))==0,img)):x_to_tip{img}(end),splined_y{img},polyfit_degree);
    pfit_y{img} = polyval(pfit,min(xf(isnan(xf(:,img))==0,img)):x_to_tip{img}(end));

    
%     % Spline the results
%     splined_y{img} = spline([xf(:,img);x_to_tip{img}],[yf(:,img);y_to_tip{img}],min(xf(isnan(xf(:,img))==0,img)):x_to_tip{img}(end));
%     splined_x{img} = min(xf(isnan(xf(:,img))==0,img)):x_to_tip{img}(end);
end

disp('All finished! Thank you!')
uicontrol('Style','text','Position',[400 200 200 75],'String','All Finished! Thanks!'); % show image number
pause(0.5);clf;
uicontrol('Style','pushbutton','Position',[400 200 200 75],...
    'String','Spline & See Results?','Callback',{@LOCAL_show_splined_results,pfit_x,pfit_y});

% title(sprintf('click Next Image to Start!'))
% h_nextimage = uicontrol('Style','pushbutton','Position',[200 300 100 50],...
%     'String','Next Image','Callback',{@LOCAL_add_points,xf,yf,prefix,fig,num_images});%{'add_points_sub',img,xf,yf })  % Control button to move to next image

function LOCAL_show_splined_results(hObject,eventdata,splined_x,splined_y)
%% function LOCAL_show_splined_results(x_to_tip,y_to_tip)
    global prefix
    figure(10201)
    for spl_img = 1:length(splined_x)
%         thisnum = sprintf('%06d',ImgNum); % pad image number
        imshow(imread([prefix,sprintf('%06d',spl_img),'.tif'])); % concatenate prefix and image number... print picture: rat
        hold on;
        uicontrol('Style','text','Position',[400 50 60 20],'String',['Frame: ',num2str(spl_img)]); % show image number
        plot(splined_x{spl_img},splined_y{spl_img});
        pause(0.05); clf
    end

function LOCAL_add_points(hObject,eventdata,xf,yf,fig,num_images)
%% function [x_to_tip,y_to_tip] = add_points_sub(stuffz)
    global ImgNum
    global x_to_tip
    global y_to_tip
    if ImgNum < num_images
       ImgNum = ImgNum + 1;
        thisnum = sprintf('%06d',ImgNum); % pad image number
        imshow(imread([prefix,thisnum,'.tif'])); % concatenate prefix and image number... print picture: rat
        title(sprintf('click on exta points\nPress ENTER to end\nPress BACKSPACE to delete previous point')) % tell the people what they need to know
        uicontrol('Style','text','Position',[400 50 60 20],'String',['Frame: ',num2str(ImgNum)]); % show image number
        hold on % stay!
        plot(xf(:,ImgNum),yf(:,ImgNum)) % plot the already plotted whisker 
        [x_to_tip{ImgNum},y_to_tip{ImgNum}] = getpts(fig);
    else
        disp('All finished! Thank you!')
        uicontrol('Style','text','Position',[400 200 200 75],'String','All Finished! Thanks!'); % show image number
    end