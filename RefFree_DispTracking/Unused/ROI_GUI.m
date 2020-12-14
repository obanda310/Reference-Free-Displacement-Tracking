function roiGUIHandle = ROI_GUI(experiment)
    % Create figure window and GUI handle roiGUIHandle
    roiGUIHandle = figure('Position',[50,100,1800,900]);
    handles = guihandles(roiGUIHandle);
    
    roiBounds = [1,1,size(experiment.images,1),size(experiment.images,2)];
    roi = experiment.images;
    setappdata(roiGUIHandle,'ROIBounds',roiBounds);
    setappdata(roiGUIHandle,'ROI',roi);
    
    noImgs = size(experiment.images,3);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Displaying raw image stack
    
    handles.rawImgAxes = axes(...
        'Units','Pixels',...
        'Position',[350,200,600,600]);
    handles.rawImgSlider = uicontrol( ...
        'Style','slider', ...
        'Min',1,'Max',noImgs, ...
        'Value',1, ...
        'SliderStep',[1/noImgs,5/noImgs], ...
        'Position',[350,150,600,20], ...
        'Callback',@rawImgSliderCallback);
    handles.rawImgSliderListener = addlistener( ...
        handles.rawImgSlider, ...
        'Value','PostSet',@rawImgSliderCallback);
    handles.rawImgSliderDisplay = uicontrol( ...
        'Style','edit', ...
        'String','1', ...
        'Position',[350,125,600,20], ...
        'Callback',@rawImgSliderCallback);
    
    function rawImgSliderCallback(~,~)
        thisFrame = round((get(handles.rawImgSlider,'Value')));
        set(handles.rawImgSliderDisplay,'String',num2str(thisFrame));
        axes(handles.rawImgAxes)
        imshow(experiment.images(:,:,thisFrame),[])

        guidata(roiGUIHandle,handles);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Selecting a region of interest
    
    handles.cropCheckbox = uicontrol( ...
        'Style','checkbox', ...
        'String','Select ROI', ...
        'Position',[20,50,150,20], ...
        'Callback',@cropCheckboxCallback);
    handles.roiImgAxes = axes(...
        'Units','Pixels',...
        'Position',[1000,200,600,600]);
    handles.roiImgSlider = uicontrol( ...
        'Style','slider', ...
        'Min',1,'Max',noImgs, ...
        'Value',1, ...
        'SliderStep',[1/noImgs,5/noImgs], ...
        'Position',[1000,150,600,20], ...
        'Callback',@roiImgSliderCallback);
    handles.roiImgSliderListener = addlistener( ...
        handles.roiImgSlider, ...
        'Value','PostSet',@roiImgSliderCallback);
    handles.roiImgSliderDisplay = uicontrol( ...
        'Style','edit', ...
        'String','1', ...
        'Position',[1000,125,600,20], ...
        'Callback',@roiImgSliderCallback);
    handles.stackStart = uicontrol( ...
        'Style','Edit', ...
        'String','1', ...
        'Position',[1000,90,60,20], ...
        'Callback',@newRangeCallback);
    handles.stackEnd = uicontrol( ...
        'Style','Edit', ...
        'String',num2str(noImgs), ...
        'Position',[1500,90,60,20], ...
        'Callback',@newRangeCallback);
    
    function newRangeCallback(src,~)
        firstFrame = str2double(get(handles.stackStart,'String'));
        lastFrame = str2double(get(handles.stackEnd,'String'));
        
        set(handles.roiImgSlider,'Min',firstFrame);
        set(handles.rawImgSlider,'Min',firstFrame);
        set(handles.roiImgSlider,'Max',lastFrame);
        set(handles.rawImgSlider,'Max',lastFrame);
        
        if src == handles.stackStart
            set(handles.roiImgSlider,'Value',firstFrame);
            set(handles.rawImgSlider,'Value',firstFrame);
        else
            set(handles.roiImgSlider,'Value',lastFrame);
            set(handles.rawImgSlider,'Value',lastFrame);
        end
        
        newNoImgs = lastFrame - firstFrame + 1;
        set(handles.rawImgSlider,'SliderStep',[1/newNoImgs,5/newNoImgs]);
        set(handles.roiImgSlider,'SliderStep',[1/newNoImgs,5/newNoImgs]);
        
        guidata(roiGUIHandle,handles);
    end
    
    function cropCheckboxCallback(~,~)
        checkStatus = get(handles.cropCheckbox,'Value');
        if checkStatus == 1 % checked
            handles.r = imrect(handles.rawImgAxes);
            handles.r.addNewPositionCallback(@roiImgSliderCallback)
        else
            delete(handles.r)
        end
        
        guidata(roiGUIHandle,handles);
    end
    function roiImgSliderCallback(~,~)
        roiBounds = round(getPosition(handles.r));
        noImgs = size(experiment.images,3);
        roi = zeros(roiBounds(4)+1,roiBounds(3)+1,noImgs);
        for i = 1:noImgs
            newImg = imcrop(experiment.images(:,:,i),roiBounds);
            roi(:,:,i) = newImg;
        end
        
        thisFrame = round((get(handles.roiImgSlider,'Value')));
        set(handles.roiImgSliderDisplay,'String',num2str(thisFrame));
        axes(handles.roiImgAxes)
        imshow(roi(:,:,thisFrame),[])

        setappdata(roiGUIHandle,'ROIBounds',roiBounds);
        setappdata(roiGUIHandle,'ROI',roi);
        guidata(roiGUIHandle,handles);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end