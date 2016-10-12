function ShowStack(InputStack)
global noImgs;
global centroids2;

hFig = figure('Position',[100 100 800 800]);

handles.axes1 = axes('Units','normalized','Position',[0 0 1 1]);

%// Create slider and listener object for smooth visualization
handles.SliderFrame = uicontrol('Units','normalized','Style','slider','Position',[.300 .050 .400 .020],'Min',1,'Max',noImgs,'Value',1,'SliderStep',[1/noImgs 2/noImgs],'Callback',@XSliderCallback);
handles.SliderxListener = addlistener(handles.SliderFrame,'Value','PostSet',@(s,e) XListenerCallBack);

handles.Text1 = uicontrol('Units','normalized','Style','Text','Position',[.400 .010 .060 .030],'String','Current frame');
handles.Edit1 = uicontrol('Units','normalized','Style','Edit','Position',[.470 .010 .100 .030],'String','1');

handles.pushbuttonPlot = uicontrol('Units','normalized','Style','pushbutton','Position',[.200 .010 .060 .030],'String','Plot Centroids','Callback',@plotCircles);

%// Use setappdata to store the image stack and in callbacks, use getappdata to retrieve it and use it. Check the docs for the calling syntax.

setappdata(hFig,'MyMatrix',InputStack); %// You could use %//setappdata(0,'MyMatrix',MyMatrix) to store in the base workspace. 

%// Display 1st frame
imshow(InputStack(:,:,1),[])

%// IMPORTANT. Update handles structure.
guidata(hFig,handles);

%// Listener callback, executed when you drag the slider.

    function XListenerCallBack

        %// Retrieve handles structure. Used to let MATLAB recognize the
        %// edit box, slider and all UI components.
        handles = guidata(gcf);

%// Here retrieve MyMatrix using getappdata.
InputStack = getappdata(hFig,'MyMatrix');

        %// Get current frame
        CurrentFrame = round((get(handles.SliderFrame,'Value')));
        set(handles.Edit1,'String',num2str(CurrentFrame));

        %// Display appropriate frame.
        imshow(InputStack(:,:,CurrentFrame),'Parent',handles.axes1);

        guidata(hFig,handles);
    end


%// Slider callback; executed when the slider is release or you press
%// the arrows.
    function XSliderCallback(~,~)

        handles = guidata(gcf);

%// Here retrieve MyMatrix using getappdata.
    InputStack = getappdata(hFig,'MyMatrix');

        CurrentFrame = round((get(handles.SliderFrame,'Value')));
        set(handles.Edit1,'String',num2str(CurrentFrame));

        imshow(InputStack(:,:,CurrentFrame),[],'Parent',handles.axes1);

        guidata(hFig,handles);
    end

    function plotCircles(~,~)
        handles = guidata(gcf);
        InputStack = getappdata(hFig,'MyMatrix');
        CurrentFrame = round((get(handles.SliderFrame,'Value')));
        set(handles.Edit1,'String',num2str(CurrentFrame));
        xs = (centroids2{CurrentFrame,1}(:,1));
        ys = (centroids2{CurrentFrame,1}(:,2));
        imshow(InputStack(:,:,CurrentFrame),[],'Parent',handles.axes1);
        hold on
        plot(xs,ys,'o')
        guidata(hFig,handles);
    end 
    


end