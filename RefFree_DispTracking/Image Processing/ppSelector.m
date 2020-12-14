function ppOptions = ppSelector()
H = 150;
h.f = figure('units','pixels','name','Select Pre-Processing Options','position',[800,480,250,H],...
    'toolbar','none','menu','none');

mTextBox = uicontrol('style','text','position',[0,H-15,200,15]);
set(mTextBox,'String','Select Pre-Processing Options');

% Create yes/no checkboxes

h.c(1) = uicontrol('style','checkbox','units','pixels',...
    'position',[10,H-35,300,15],'string','Use Remove Large v2','Value',0);
h.c(2) = uicontrol('style','checkbox','units','pixels',...
    'position',[10,H-35,300,15],'string','Use Remove Large v2','Value',0);


h.diameterLabel = uicontrol('Style','Text','Units','pixels','Position',[10,H-55,300,15], ...
    'String','Approximate Feature XY Width (pixels):','HorizontalAlignment','Left');

h.diameterInput = uicontrol( ...
    'Style','Edit', ...
    'String','7.0', ...
    'Units','pixels', ...
    'Position',[75,H-75,100,15], ...
    'Callback',@selectDiameter);

newDxy = 7.0;

h.thresholdLabel = uicontrol( ...
    'Style','Text', ...
    'Units','pixels','HorizontalAlignment','Left', ...
    'Position',[10,H-95,300,15], ...
    'String','Approximate Feature Z Length(pixels):');

h.thresholdInput = uicontrol( ...
    'Style','Edit', ...
    'String','11', ...
    'Units','pixels', ...
    'Position',[75,H-115,100,15], ...
    'Callback',@selectThreshold);

newDz = 11;




    function selectDiameter(~,~)
        % Get value from h.diameterInput textbox
        newDxy = str2double(get(h.diameterInput,'String'));
    end

    function selectThreshold(~,~)
        % Get value from h.thresholdInput textbox
        newDz = str2double(get(h.thresholdInput,'String'));
    end


% Create OK pushbutton
h.p = uicontrol('style','pushbutton','units','pixels',...
    'position',[90,H-140,70,20],'string','OK',...
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
disp(newDxy)
disp(newDz)
ppOptions{1} = checked;
ppOptions{2} = newDxy;
ppOptions{3} = newDz;
end