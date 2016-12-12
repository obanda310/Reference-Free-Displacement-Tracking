function ppOptions = ppSelector()
h.f = figure('units','pixels','name','Select Pre-Processing Options','position',[800,480,300,205],...
    'toolbar','none','menu','none');

mTextBox = uicontrol('style','text','position',[0,190,200,15]);
set(mTextBox,'String','Select Pre-Processing Options');

% Create yes/no checkboxes
h.c(1) = uicontrol('style','checkbox','units','pixels',...
    'position',[10,30,300,15],'string','Subtract 95%ile Final Frame Value');
h.c(2) = uicontrol('style','checkbox','units','pixels',...
    'position',[10,50,300,15],'string','Use DoG Filter (option not yet implemented)');

h.diameterLabel = uicontrol( ...
    'Style','Text', ...
    'Units','pixels', ...
    'Position',[10,130,300,15], ...
    'String','Enter Approximate Feature Diameter for DoG (Microns)');

h.diameterInput = uicontrol( ...
    'Style','Edit', ...
    'String','1.2', ...
    'Units','pixels', ...
    'Position',[100,110,100,15], ...
    'Callback',@selectDiameter);

newD = 1.2;

h.thresholdLabel = uicontrol( ...
    'Style','Text', ...
    'Units','pixels', ...
    'Position',[10,90,300,15], ...
    'String','Enter Feature Threshold % After DoG');

h.thresholdInput = uicontrol( ...
    'Style','Edit', ...
    'String','25', ...
    'Units','pixels', ...
    'Position',[100,70,100,15], ...
    'Callback',@selectThreshold);

newT = 25;

    function selectDiameter(~,~)
% Get value from h.diameterInput textbox
        newD = str2double(get(h.diameterInput,'String'));
    end

    function selectThreshold(~,~)
% Get value from h.thresholdInput textbox
        newT = str2double(get(h.thresholdInput,'String'));
   end

% Create OK pushbutton
h.p = uicontrol('style','pushbutton','units','pixels',...
    'position',[40,5,70,20],'string','OK',...
    'callback',@p_call);
uiwait(h.f);


% Pushbutton callback
    function p_call(varargin)
        vals = get(h.c,'Value');
        checked = find([vals{:}]);
        if isempty(checked)
            checked = 'none';
        end
        close;
    end
disp(checked)
disp(newD)
ppOptions{1} = checked;
ppOptions{2} = newD;
ppOptions{3} = newT;
end