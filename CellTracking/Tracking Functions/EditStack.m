% InputStack should be a 3D image matrix where dimensions 1 and 2 contain
% pixel information in y and x, respectively, and dimension 3 indexes
% separate images in the stack
function [OutputStackFinal,cropImages,bounds] = EditStack(InputStack,Original,centroids)
    
    noImgs = size(InputStack,3);
    bounds = [1,1,size(InputStack,1),size(InputStack,2)];
    OutputStack = InputStack;
    cropImages = Original;
    hFig = figure('Position',[100 100 800 800]);

    handles.axes1 = axes('Units','normalized','Position',[0 0 1 1]);

    % Create slider and listener object for smooth visualization
    handles.SliderFrame = uicontrol( ...
        'Style','slider', ...
        'Units','normalized', ... 
        'Position',[.300 .050 .400 .020], ...
        'Min',1, ...
        'Max',noImgs, ...
        'Value',1, ...
        'SliderStep',[1/noImgs 2/noImgs], ...
        'Callback',@XSliderCallback);
    handles.SliderxListener = addlistener( ...
        handles.SliderFrame,'Value','PostSet',@(s,e) XListenerCallBack);

    handles.Text1 = uicontrol( ...
        'Style','Text', ...
        'Units','normalized', ...
        'Position',[.400 .010 .060 .030], ...
        'String','Current frame');
    handles.Edit1 = uicontrol( ...
        'Style','Edit', ...
        'Units','normalized', ...
        'Position',[.470 .010 .100 .030], ...
        'String','1');
    
    % Plot Centroids2 over current image
    handles.pushbuttonPlot = uicontrol( ...
        'Style','pushbutton', ...
        'Units','normalized', ...
        'Position',[.100 .010 .130 .030], ...
        'String','Plot Centroids', ...
        'Callback',@plotCircles);
    
    % Crop current image and update stack information    
    handles.pushbuttonPlot = uicontrol( ...
        'Style','pushbutton', ...
        'Units','normalized', ...
        'Position',[.100 .045 .060 .030], ...
        'String','Crop', ...
        'Callback',@cropImage);
    
    % Undo Crop
    handles.pushbuttonPlot = uicontrol( ...
        'Style','pushbutton', ...
        'Units','normalized', ...
        'Position',[.165 .045 .060 .030], ...
        'String','Original', ...
        'Callback',@revert);
    
    handles.pushbuttonPlot = uicontrol( ...
        'Style','pushbutton', ...
        'Units','normalized', ...
        'Position',[.1 .08 .13 .030], ...
        'String','Accept and Close', ...
        'Callback',@acceptandclose);
    

    % Use setappdata to store the image stack and in callbacks, use 
    % getappdata to retrieve it and use it. Check the docs for the calling 
    % syntax.
    % You could use setappdata(0,'MyMatrix',MyMatrix) to store in the base 
    % workspace. 
    setappdata(hFig,'MyMatrix',OutputStack);
    setappdata(hFig,'MyOrigMatrix',cropImages);

    % Display 1st frame
    imshow(OutputStack(:,:,1),[])

    % IMPORTANT. Update handles structure.
    guidata(hFig,handles);
    waitfor(hFig)

    % Listener callback, executed when you drag the slider.
    function XListenerCallBack
        % Retrieve handles structure. Used to let MATLAB recognize the
        % edit box, slider and all UI components.
        handles = guidata(gcf);
        % Here retrieve MyMatrix using getappdata.
        OutputStack = getappdata(hFig,'MyMatrix');
        % Get current frame
        CurrentFrame = round((get(handles.SliderFrame,'Value')));
        set(handles.Edit1,'String',num2str(CurrentFrame));
        % Display appropriate frame.
        imshow(OutputStack(:,:,CurrentFrame),[]);
        guidata(hFig,handles);
    end

    % Slider callback; executed when the slider is release or you press
    % the arrows.
    function XSliderCallback(~,~)
        handles = guidata(gcf);
        % Here retrieve MyMatrix using getappdata.
        OutputStack = getappdata(hFig,'MyMatrix');
        CurrentFrame = round((get(handles.SliderFrame,'Value')));
        set(handles.Edit1,'String',num2str(CurrentFrame));
        imshow(OutputStack(:,:,CurrentFrame),[]);
        guidata(hFig,handles);
    end

    function plotCircles(~,~)
        handles = guidata(gcf);
        OutputStack = getappdata(hFig,'MyMatrix');
        CurrentFrame = round((get(handles.SliderFrame,'Value')));
        set(handles.Edit1,'String',num2str(CurrentFrame));
        xs = (centroids{CurrentFrame,1}(:,1));
        ys = (centroids{CurrentFrame,1}(:,2));
        imshow(OutputStack(:,:,CurrentFrame),[]);
        hold on
        plot(xs,ys,'o')
        guidata(hFig,handles);
        hold off
    end 
    
    function cropImage(~,~)
        f = msgbox('Double-Click to crop after making rectangular selection!');
        waitfor(f);
        handles = guidata(gcf);
        OutputStack = getappdata(hFig,'MyMatrix');
        CurrentFrame = round((get(handles.SliderFrame,'Value')));
        set(handles.Edit1,'String',num2str(CurrentFrame));    
        
        [~, bounds] = imcrop(OutputStack(:,:,CurrentFrame));
        bounds = round(bounds);
        CroppedStack = OutputStack(bounds(1,2):bounds(1,2)+bounds(1,4),bounds(1,1):bounds(1,1)+bounds(1,3),:);
        cropImages = Original(bounds(1,2):bounds(1,2)+bounds(1,4),bounds(1,1):bounds(1,1)+bounds(1,3),:);
        setappdata(hFig,'MyMatrix',CroppedStack);
        setappdata(hFig,'MyOrigMatrix',cropImages);
        OutputStack = getappdata(hFig,'MyMatrix');
        imshow(OutputStack(:,:,CurrentFrame),[],'Parent',handles.axes1);
        guidata(hFig,handles);
    end

    function revert(~,~)
        handles = guidata(gcf);
        OutputStack = getappdata(hFig,'MyMatrix');
        CurrentFrame = round((get(handles.SliderFrame,'Value')));
        set(handles.Edit1,'String',num2str(CurrentFrame));
        setappdata(hFig,'MyMatrix',InputStack);
        OutputStack = getappdata(hFig,'MyMatrix');
        imshow(OutputStack(:,:,CurrentFrame),[]);
        guidata(hFig,handles);
    end
    
    function acceptandclose(~,~)
        OutputStackFinal = OutputStack;
        close()
    end
end