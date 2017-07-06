%GUI to be run along with trajectories.m to manually select output image
%files.
function y1 = OutputSelector()
% Create figure

h.f = figure('units','pixels','name','Select Outputs','position',[300,480,750,205],...
    'toolbar','none','menu','none');



% Create yes/no checkboxes
%%%%%%%%%%%%%%%%%%%%%%%%
%Options for plotting
%%%%%%%%%%%%%%%%%%%%%%%%

mTextBox = uicontrol('style','text','position',[0,190,250,15]);
set(mTextBox,'String','Plotting Displacement Vectors and Centroids')

h.c(2) = uicontrol('style','checkbox','units','pixels',...
    'position',[10,110,300,15],'string','Zero-State Displacement Fields on Black Background');
h.c(3) = uicontrol('style','checkbox','units','pixels',...
    'position',[10,150,300,15],'string','DEBUGGING:Plotting Corrected X,Y Coordinate Fields','Value',1);
h.c(7) = uicontrol('style','checkbox','units','pixels',...
    'position',[20,130,300,15],'string','No Centroids');
h.c(8) = uicontrol('style','checkbox','units','pixels',...
    'position',[110,130,300,15],'string','No Quiver');
h.c(4) = uicontrol('style','checkbox','units','pixels',...
    'position',[10,170,300,15],'string','Centroids on Black Background');
h.c(5) = uicontrol('style','checkbox','units','pixels',...
    'position',[10,90,300,15],'string','Transmitted/Trajectory Quiver Overlays with Save','Value',1);

%Unused (old)
mTextBox = uicontrol('style','text','position',[0,70,150,15]);
set(mTextBox,'String','Vector/Centroids Settings:')

h.c(1) = uicontrol('style','checkbox','units','pixels',...
    'position',[10,30,300,15],'string','Use Default Colormap','Value',1);
h.c(6) = uicontrol('style','checkbox','units','pixels',...
    'position',[10,50,300,15],'string','Use "Native" Image Size');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Options for Intensity value operations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mTextBox = uicontrol('style','text','position',[350,190,200,15]);
set(mTextBox,'String','Fitting Intensity Values of Pillars')
h.c(17) = uicontrol('style','checkbox','units','pixels',...
    'position',[360,170,300,15],'string','Attempt Intensity Z');
h.c(9) = uicontrol('style','checkbox','units','pixels',...
    'position',[370,130,300,15],'string','Pillar Analysis');
h.c(10) = uicontrol('style','checkbox','units','pixels',...
    'position',[370,150,300,15],'string','PlaneFit');
h.c(14) = uicontrol('style','checkbox','units','pixels',...
    'position',[370,110,300,15],'string','Profile Fits');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Options for Creating FE Meshes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mTextBox = uicontrol('style','text','position',[550,190,150,15]);
set(mTextBox,'String','Creating Mesh for FE')

h.c(11) = uicontrol('style','checkbox','units','pixels',...
    'position',[560,170,300,15],'string','Shear Mesh');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Options for Creating Heat Maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mTextBox = uicontrol('style','text','position',[290,95,200,15]);
set(mTextBox,'String','Heat Maps')

h.c(15) = uicontrol('style','checkbox','units','pixels',...
    'position',[360,80,300,15],'string','XY Heatmap','Value',1);
h.c(13) = uicontrol('style','checkbox','units','pixels',...
    'position',[370,60,300,15],'string','Through-Z Normal');

h.c(12) = uicontrol('style','checkbox','units','pixels',...
    'position',[370,60,300,15],'string','Through-Z Normal');

h.c(16) = uicontrol('style','checkbox','units','pixels',...
    'position',[360,20,300,15],'string','Z Heatmap','Value',1);



% Create OK pushbutton
h.p = uicontrol('style','pushbutton','units','pixels',...
    'position',[40,5,70,20],'string','OK',...
    'callback',@p_call);
uicontrol(h.p)
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