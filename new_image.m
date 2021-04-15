function new_image(hObject,eventdata,xf,yf,prefix,fig,num_images)
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