%GUI to be run along with trajectories.m to manually select output image
%files.
function y1 = OutputSelector()
% Create figure

h.f = figure('units','pixels','name','Select Outputs','position',[800,480,300,205],...
    'toolbar','none','menu','none');

mTextBox = uicontrol('style','text','position',[0,190,200,15])
set(mTextBox,'String','Fitting Intensity Values of Pillars')
mTextBox = uicontrol('style','text','position',[0,130,250,15])
set(mTextBox,'String','Plotting Displacement Vectors and Centroids')

% Create yes/no checkboxes
h.c(1) = uicontrol('style','checkbox','units','pixels',...
    'position',[10,30,300,15],'string','Zero-State Displacement Fields, Old Version');
h.c(2) = uicontrol('style','checkbox','units','pixels',...
    'position',[10,50,300,15],'string','Zero-State Displacement Fields on Black Background');
h.c(3) = uicontrol('style','checkbox','units','pixels',...
    'position',[10,90,300,15],'string','DEBUGGING:Plotting Corrected X,Y Coordinate Fields');
h.c(7) = uicontrol('style','checkbox','units','pixels',...
    'position',[20,70,300,15],'string','No Centroids');
h.c(8) = uicontrol('style','checkbox','units','pixels',...
    'position',[110,70,300,15],'string','No Quiver');
h.c(4) = uicontrol('style','checkbox','units','pixels',...
    'position',[10,110,300,15],'string','Centroids on Black Background');
h.c(5) = uicontrol('style','checkbox','units','pixels',...
    'position',[10,30,300,15],'string','Transmitted/Trajectory Quiver Overlays with Save');
h.c(6) = uicontrol('style','checkbox','units','pixels',...
    'position',[10,150,300,15],'string','Transmitted/Trajectory Quiver Overlays without Save');
h.c(9) = uicontrol('style','checkbox','units','pixels',...
    'position',[10,150,300,15],'string','PillarFit');
h.c(10) = uicontrol('style','checkbox','units','pixels',...
    'position',[10,170,300,15],'string','PlaneFit');
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
y1 = checked;
end