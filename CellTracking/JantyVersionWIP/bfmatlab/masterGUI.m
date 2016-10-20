function masterHandle = masterGUI
    % Create figure window and GUI handle h
    masterHandle = figure('Position',[50,100,1800,900]);
    
    handles = guihandles(masterHandle);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Choosing .czi file to be analyzed
    
    handles.chooseCziButton = uicontrol( ...
        'Style','pushbutton', ...
        'String','Choose .czi file', ...
        'Position',[20,20,150,20], ...
        'Callback',@chooseCziButtonCallback);
    
    function chooseCziButtonCallback(~,~)
        [images,meta] = getCzi;
        
        axes(handles.rawImgAxes)
        imshow(images(:,:,1),[])
        
        noImgs = size(images,3);
        set(handles.rawImgSlider, ...
            'Max',noImgs, ...
            'SliderStep',[1/noImgs,5/noImgs]);
        
        setappdata(masterHandle,'images',images);
        setappdata(masterHandle,'metadata',meta);
        guidata(masterHandle,handles);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Displaying raw image stack
    
    handles.rawImgAxes = axes(...
        'Units','Pixels',...
        'Position',[350,200,600,600]);
    handles.rawImgSlider = uicontrol( ...
        'Style','slider', ...
        'Min',1,'Max',50, ...
        'Value',1, ...
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
        images = getappdata(masterHandle,'images');
        thisFrame = round((get(handles.rawImgSlider,'Value')));
        set(handles.rawImgSliderDisplay,'String',num2str(thisFrame));
        axes(handles.rawImgAxes)
        imshow(images(:,:,thisFrame),[])

        guidata(masterHandle,handles);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Displaying raw image stack
end