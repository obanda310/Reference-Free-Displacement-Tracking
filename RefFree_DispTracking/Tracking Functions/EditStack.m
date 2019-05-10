% InputStack should be a 3D image matrix where dimensions 1 and 2 contain
% pixel information in y and x, respectively, and dimension 3 indexes
% separate images in the stack
function [roiMaskStack,roiImgStack,roiBounds,redoCheck] = EditStack(maskStack,originalStack,expCheck) %,centroids

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Initialize GUI

% Create GUI figure window
hFig = figure('Position',[100 100 800 800]);

iLow = 0;
iHigh = prctile(maskStack(:),95)*2;
contrast = [iLow iHigh];

% noImgs is the number of images in the input stack
noImgs = size(maskStack,3);
% Initialize the bounds of the ROI to be equal to the size of an entire
% image
roiBounds = [1,1,size(maskStack,2),size(maskStack,1)];

firstFrame = 1;
lastFrame = noImgs;
redoCheck = 0;
% Initialize the OutputStack to be equal to the input mask stack
OutputStack = maskStack;
% Initialize the ROI image stack to be equal to the original image
% stack
roiImgStack = originalStack;

% Use setappdata to store the image stack and in callbacks, use
% getappdata to retrieve it and use it. Check the docs for the calling
% syntax.
% You could use setappdata(0,'MyMatrix',MyMatrix) to store in the base
% workspace.
setappdata(hFig,'MyMatrix',OutputStack);
setappdata(hFig,'MyOrigMatrix',roiImgStack);

% Create axes for the images to be displayed in
handles.axes1 = axes( ...
    'Units','normalized', ...
    'Position',[0 0 1 1]);

% Display 1st frame
imshow(OutputStack(:,:,1),[contrast])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Create slider and callback to scroll through image stack

% Create slider and listener object for smooth visualization
handles.SliderFrame = uicontrol( ...
    'Style','slider', ...
    'Units','normalized', ...
    'Position',[.300 .060 .400 .020], ...
    'Min',1, ...
    'Max',noImgs, ...
    'Value',1, ...
    'SliderStep',[1/noImgs 2/noImgs], ...
    'Callback',@XSliderCallback);
% Slider listener calls the XListenerCallback immediately after the
% slider position is changed - this is what allows smooth visualization
handles.SliderxListener = addlistener( ...
    handles.SliderFrame,'Value','PostSet',@(s,e) XListenerCallback);

% Slider callback; executed when the slider is release or you press
% the arrows.
    function XSliderCallback(~,~)
        handles = guidata(gcf);
        % Here retrieve MyMatrix using getappdata.
        OutputStack = getappdata(hFig,'MyMatrix');
        CurrentFrame = round((get(handles.SliderFrame,'Value')));
        set(handles.currentFrameSelect,'String',num2str(CurrentFrame));
        imshow(OutputStack(:,:,CurrentFrame),[contrast]);
        guidata(hFig,handles);
    end

% Listener callback, executed when you drag the slider.
    function XListenerCallback
        % Retrieve handles structure. Used to let MATLAB recognize the
        % edit box, slider and all UI components.
        handles = guidata(gcf);
        % Here retrieve MyMatrix using getappdata.
        OutputStack = getappdata(hFig,'MyMatrix');
        % Get current frame
        CurrentFrame = round((get(handles.SliderFrame,'Value')));
        set(handles.currentFrameSelect,'String',num2str(CurrentFrame));
        % Display appropriate frame.
        imshow(OutputStack(:,:,CurrentFrame),[contrast]);
        guidata(hFig,handles);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Allow use to select start and end images in stack

% Create label for textbox in which user can select the image number
% they would like to be the start of the image stack
handles.stackStartLabel = uicontrol( ...
    'Style','Text', ...
    'Units','normalized', ...
    'Position',[.300 .035 .100 .020], ...
    'String','First Frame');
% Create textbox in which user can select the image number they would
% like to be the start of the image stack
handles.stackStart = uicontrol( ...
    'Style','Edit', ...
    'String','1', ...
    'Units','normalized', ...
    'Position',[.300,.010,.100,.020], ...
    'Callback',@newRangeCallback);
% Create label for textbox in which user can select the image number
% they would like to be the end of the image stack
handles.stackEndLabel = uicontrol( ...
    'Style','Text', ...
    'Units','normalized', ...
    'Position',[.600 .035 .100 .020], ...
    'String','Last Frame');
% Create textbox in which user can select the image number they would
% like to be the end of the image stack
handles.stackEnd = uicontrol( ...
    'Style','Edit', ...
    'String',num2str(noImgs), ...
    'Units','normalized', ...
    'Position',[.600,.010,.100,.020], ...
    'Callback',@newRangeCallback);

% Callback will adjust slider range based upon the user-input values
% for start frame and end frame. Will also immediately display the
% image corresponding to the frame the user just specified.
    function newRangeCallback(src,~)
        % Get value from handles.stackStart textbox
        firstFrame = str2double(get(handles.stackStart,'String'));
        % Get value from handles.stackEnd textbox
        lastFrame = str2double(get(handles.stackEnd,'String'));
        
        % Set the minimum slider value to be equal to the value obtained
        % from handles.stackStart
        set(handles.SliderFrame,'Min',firstFrame);
        % Set the maximum slider value to be equal to the value obtained
        % from handles.stackEnd
        set(handles.SliderFrame,'Max',lastFrame);
        
        % If the user just edited the value in the handles.stackStart
        % textbox, set the slider to that frame
        if src == handles.stackStart
            set(handles.SliderFrame,'Value',firstFrame);
            % Otherwise (i.e. if the user just edited the value in the
            % handles.stackEnd textbox), set the slider to that frame
        else
            set(handles.SliderFrame,'Value',lastFrame);
        end
        % The above if-statement is necessary because when changing the
        % slider bounds, if the slide position is set to a value outside of
        % those bounds, there will be an error. Therefore, we move the
        % slider to be within the bounds of the image by setting the slider
        % to the frame the user just selected as a minimum or maximum
        % bound for the image stack
        
        % Determine the new number of images in the stack based upon the
        % first frame and last frame set by the user, and adjust the slider
        % step-size accordingly
        newNoImgs = lastFrame - firstFrame + 1;
        set(handles.SliderFrame,'SliderStep',[1/newNoImgs,2/newNoImgs]);
        
        setappdata(hFig,'FirstFrame',firstFrame)
        setappdata(hFig,'LastFrame',lastFrame)
        guidata(hFig,handles);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Show current frame number and allow user to set manually

% Create label for textbox in which user can set the image number they
% would like to view
handles.currentFrameLabel = uicontrol( ...
    'Style','Text', ...
    'Units','normalized', ...
    'Position',[.450 .035 .100 .020], ...
    'String','Current Frame');
% Create textbox in which user can set the image number they would like
% to view
handles.currentFrameSelect = uicontrol( ...
    'Style','Edit', ...
    'Units','normalized', ...
    'Position',[.450 .010 .100 .020], ...
    'String','1', ...
    'Callback',@currentFrameSelectCallback);

% Callback will change the slider position to the image number the user
% specifies in the handles.currentFrameSelect textbox. Upon moving the
% slider, the slider callback functions will update the image being
% displayed
    function currentFrameSelectCallback(~,~)
        % Get value of user-input frame from the handles.currentFrameSelect
        % textbox
        selectedFrame = str2double(get(handles.currentFrameSelect,'String'));
        % Set the slider position to the value in the
        % handles.currentFrameSelect textbox
        set(handles.SliderFrame,'Value',selectedFrame);
        
        guidata(hFig,handles);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Allow user to select ROI sub-section of image stack

%%%%%%%%%%%%%%%%%%
%%%% CROPPING %%%%
%%%%%%%%%%%%%%%%%%

% Create button that, when clicked, will use cropImgCallback to allow
% the user to select an ROI
handles.cropImgButton = uicontrol( ...
    'Style','pushbutton', ...
    'Units','normalized', ...
    'Position',[.100 .035 .060 .020], ...
    'String','Crop', ...
    'Callback',@cropImgCallback);

% Callback will open a tool to create a resizeable and movable
% rectangle for the user to select an ROI
    function cropImgCallback(~,~)
        % Instruct user how to select ROI
        f = msgbox('Double-Click to crop after making rectangular selection!');
        % Wait for user to click "OK" before allowing the user to proceed
        % to select ROI
        waitfor(f);
        
        OutputStack = getappdata(hFig,'MyMatrix');
        CurrentFrame = round((get(handles.SliderFrame,'Value')));
        
        % imcrop allows user to select a rectangular ROI and outputs the
        % bounds of the ROI as four-element vector with the following
        % values: [xmin,ymin,width,height]
        [~, roiBounds] = imcrop;
        roiBounds = round(roiBounds);
        % Minimum x-position of ROI bounding rectangle; note, x-position
        % corresponds to column index
        roiX = roiBounds(1);
        % Minimum y-position of ROI bounding rectangle; note, y-position
        % corresponds to row index
        roiY = roiBounds(2);
        % Width of ROI bounding rectangle (spanning columns)
        roiW = roiBounds(3);
        % Height of ROI bounding rectangle (spanning rows)
        roiH = roiBounds(4);
        
        CroppedStack = OutputStack(roiY:roiY+roiH,roiX:roiX+roiW,:);
        roiImgStack = originalStack(roiBounds(2):roiBounds(2)+roiBounds(4),roiBounds(1):roiBounds(1)+roiBounds(3),:);
        setappdata(hFig,'MyMatrix',CroppedStack);
        setappdata(hFig,'MyOrigMatrix',roiImgStack);
        
        OutputStack = getappdata(hFig,'MyMatrix');
        imshow(OutputStack(:,:,CurrentFrame),[],'Parent',handles.axes1);
        guidata(hFig,handles);
    end

%%%%%%%%%%%%%%%%%%%%%%%
%%%% UNDO CROPPING %%%%
%%%%%%%%%%%%%%%%%%%%%%%

% Undo Crop
handles.cropImgButton = uicontrol( ...
    'Style','pushbutton', ...
    'Units','normalized', ...
    'Position',[.165 .035 .060 .020], ...
    'String','Original', ...
    'Callback',@revert);

    function revert(~,~)
        handles = guidata(gcf);
        OutputStack = getappdata(hFig,'MyMatrix');
        CurrentFrame = round((get(handles.SliderFrame,'Value')));
        set(handles.currentFrameSelect,'String',num2str(CurrentFrame));
        setappdata(hFig,'MyMatrix',maskStack);
        OutputStack = getappdata(hFig,'MyMatrix');
        imshow(OutputStack(:,:,CurrentFrame),[contrast]);
        guidata(hFig,handles);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% FINALIZE CROPPING %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

handles.cropImgButton = uicontrol( ...
    'Style','pushbutton', ...
    'Units','normalized', ...
    'Position',[.100 .060 .125 .020], ...
    'String','Accept and Close', ...
    'Callback',@acceptandclose);

    function acceptandclose(~,~)
        roiMaskStack = OutputStack(:,:,firstFrame:lastFrame);
        roiImgStack = roiImgStack(:,:,firstFrame:lastFrame);
        close()
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Redo Preprocessing with New Parameters %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if expCheck == 1
    handles.redoPP = uicontrol( ...
        'Style','pushbutton', ...
        'Units','normalized', ...
        'Position',[.100 .00 .125 .020], ...
        'String','Redo Preprocessing', ...
        'Callback',@redoPP);
end

    function redoPP(~,~)
        roiMaskStack = OutputStack(:,:,firstFrame:lastFrame);
        roiImgStack = roiImgStack(:,:,firstFrame:lastFrame);
        redoCheck = 1;
        close()
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%% Plot 2D centroids
%
%     % Plot Centroids2 over current image
%     handles.cropImgButton = uicontrol( ...
%         'Style','pushbutton', ...
%         'Units','normalized', ...
%         'Position',[.100 .010 .130 .030], ...
%         'String','Plot Centroids', ...
%         'Callback',@plotCircles);
%
%     function plotCircles(~,~)
%         handles = guidata(gcf);
%         OutputStack = getappdata(hFig,'MyMatrix');
%         CurrentFrame = round((get(handles.SliderFrame,'Value')));
%         set(handles.currentFrameSelect,'String',num2str(CurrentFrame));
%         xs = (centroids{CurrentFrame,1}(:,1));
%         ys = (centroids{CurrentFrame,1}(:,2));
%         imshow(OutputStack(:,:,CurrentFrame),[]);
%         hold on
%         plot(xs,ys,'o')
%         guidata(hFig,handles);
%         hold off
%     end
%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% IMPORTANT. Update handles structure.
guidata(hFig,handles);
waitfor(hFig)
end