%% Automated Script
clear all
close all
[ListName,ListPath] = uigetfile;
fullListName = strcat(ListPath,ListName);
load(strcat(ListPath,ListName),'dirList')
parts = strsplit(ListName, 'List');
prefix = parts{1};
set(0,'defaultfigurecolor',[1 1 1])
%% For updating paths (when files are moved)
%set a as new path root
%set b as old path root
ak = strfind(fullListName,ListName);

a = [fullListName(1:ak-1) ''];

Schar = char(dirList{:});
all_rows_same = all(diff(Schar, [], 1) == 0, 1);
common_cols = find(~all_rows_same, 1, 'first');
if isempty(common_cols)
    b = '?'
else
    b = dirList{1}(1:common_cols-1);
end

for i = 1:size(dirList,1)
    
    dirList{i,1} = strrep(dirList{i,1},b,a);
end
%%
save(fullListName,'dirList')

%% For clearing old files
for j = 1:2
    for i = 1:size(dirList,1)
        cleanOld(dirList{i,1})
        disp(num2str(i))
    end
end
%%
clear AFMResults
failed = 0;
for i = 1:size(dirList,1)
    try
        disp(['Trying ' num2str(i)])
        AFMResults{i} = AFMCurveRead_Plot(dirList{i,1});
    catch
        disp(['AFM Results failed at ' num2str(i)])
        failed(i) = i;
    end
end
%% Figure Settings

%CHANGE FONT SIZES HERE
AxisFontSize = 24;
AxisTitleFontSize = 16;
LegendFontSize = 16;

colOptions{1,1} = 'white';
colOptions{2,1} = 'black';
colOptions{1,2} = 'black';
colOptions{2,2} = 'white';

%% Overlay all force-indentation curves
pattern = [1 2 3];


fcolor = colOptions{1,2};
bcolor = colOptions{2,2};
PA = figure;
hold on

pTotal = zeros(size(AFMResults{1}{1},1),4);
pNum = 0;
npTotal = zeros(size(AFMResults{1}{1},1),4);
npNum = 0;
for j = nopattern
    pNum = pNum + size(AFMResults{j},2);
    for i = 1:size(AFMResults{j},2)
        pTotal(:,1) = pTotal(:,1)+AFMResults{j}{i}.Height_Sensor_nm_Ex/1000;
        pTotal(:,2) = pTotal(:,2)+AFMResults{j}{i}.Height_Sensor_nm_Rt/1000;
        pTotal(:,3) = pTotal(:,3)+AFMResults{j}{i}.Defl_pN_Ex/1000;
        pTotal(:,4) = pTotal(:,4)+AFMResults{j}{i}.Defl_pN_Rt/1000;
        plot(AFMResults{j}{i}.Height_Sensor_nm_Ex/1000,AFMResults{j}{i}.Defl_pN_Ex/1000,'color',[.7 .7 1],'linestyle','-.','linewidth',1)
        plot(AFMResults{j}{i}.Height_Sensor_nm_Rt/1000,AFMResults{j}{i}.Defl_pN_Rt/1000,'color',[1 .7 .7],'linestyle','-.','linewidth',1)
    end
end

pMean = pTotal/pNum;
plot(pMean(:,1),pMean(:,3),'black','linewidth',1)
plot(pMean(:,2),pMean(:,4),'black','linewidth',1)

% Figure Settings
ylim([-10 80])

set(gca,'Color',bcolor)
%Axes, Text, Legends
set(gca,'fontsize',AxisFontSize,'XColor',fcolor,'YColor',fcolor,'YMinorTick','on','LineWidth',2)
ytickformat('%.1f')
xt = 'Indentation (\mum)';% input('enter the xaxis label','s');
yt = 'Force (nN)'; %input('enter the yaxis label','s');
tt = 'Line-Profile Displacements';%input('enter the title','s');
%le = 'Single'; %input('enter the legend','s');
le2 = 'Time-Lapse 1';
le3 = 'Time-Lapse 2';
xl = xlabel(xt);
yl = ylabel(yt);
%tl = title(tt);

set(xl, 'fontweight','bold','fontsize',28,'color',fcolor);
set(yl,'fontweight','bold','fontsize',28,'color',fcolor);
%leg = legend(['\color{' fcolor '}' le2],['\color{' fcolor '}' le3],'location','west','fontcolor',fcolor);
%legend boxoff
%leg.FontSize = LegendFontSize;
%set(tl,'fontweight','bold','fontsize',title_font_size)
axis([-7 0 -10 80])
xticks([-6 -4 -2 0])

title = ['\Patterned ForceInd Curves ' fcolor ' on ' bcolor];
savefile = [ListPath title];
export_fig(PA,savefile,'-native');
%%
nopattern = [4 5 6];

NP = figure;
hold on
for j = pattern
    npNum = npNum + size(AFMResults{j},2);
    for i = 1:size(AFMResults{j},2)
        npTotal(:,1) = npTotal(:,1)+AFMResults{j}{i}.Height_Sensor_nm_Ex/1000;
        npTotal(:,2) = npTotal(:,2)+AFMResults{j}{i}.Height_Sensor_nm_Rt/1000;
        npTotal(:,3) = npTotal(:,3)+AFMResults{j}{i}.Defl_pN_Ex/1000;
        npTotal(:,4) = npTotal(:,4)+AFMResults{j}{i}.Defl_pN_Rt/1000;
        plot(AFMResults{j}{i}.Height_Sensor_nm_Ex/1000,AFMResults{j}{i}.Defl_pN_Ex/1000,'color',[.7 .7 1],'linestyle','-.','linewidth',1)
        plot(AFMResults{j}{i}.Height_Sensor_nm_Rt/1000,AFMResults{j}{i}.Defl_pN_Rt/1000,'color',[1 .7 .7],'linestyle','-.','linewidth',1)
    end
end

npMean = npTotal/npNum;
plot(npMean(:,1),npMean(:,3),'black','linewidth',1)
plot(npMean(:,2),npMean(:,4),'black','linewidth',1)

% Figure Settings
ylim([-10 80])

set(gca,'Color',bcolor)
%Axes, Text, Legends
set(gca,'fontsize',AxisFontSize,'XColor',fcolor,'YColor',fcolor,'YMinorTick','on','LineWidth',2)
ytickformat('%.1f')
xt = 'Indentation (\mum)';% input('enter the xaxis label','s');
yt = 'Force (nN)'; %input('enter the yaxis label','s');
tt = 'Line-Profile Displacements';%input('enter the title','s');
%le = 'Single'; %input('enter the legend','s');
le2 = 'Time-Lapse 1';
le3 = 'Time-Lapse 2';
xl = xlabel(xt);
yl = ylabel(yt);
%tl = title(tt);

set(xl, 'fontweight','bold','fontsize',28,'color',fcolor);
set(yl,'fontweight','bold','fontsize',28,'color',fcolor);
%leg = legend(['\color{' fcolor '}' le2],['\color{' fcolor '}' le3],'location','west','fontcolor',fcolor);
%legend boxoff
%leg.FontSize = LegendFontSize;
%set(tl,'fontweight','bold','fontsize',title_font_size)
axis([-7 0 -10 80])
xticks([-6 -4 -2 0])

title = ['\Non-Patterned ForceInd Curves ' fcolor ' on ' bcolor];
savefile = [ListPath title];
export_fig(NP,savefile,'-native');
%%
PaVsNp = figure;
hold on
plot(pMean(:,1),pMean(:,3),'red')
plot(pMean(:,2),pMean(:,4),'red')
plot(npMean(:,1),npMean(:,3),'blue')
plot(npMean(:,2),npMean(:,4),'blue')

% Figure Settings
ylim([-10 80])

set(gca,'Color',bcolor)
%Axes, Text, Legends
set(gca,'fontsize',AxisFontSize,'XColor',fcolor,'YColor',fcolor,'YMinorTick','on','LineWidth',2)
ytickformat('%.1f')
xt = 'Indentation (\mum)';% input('enter the xaxis label','s');
yt = 'Force (nN)'; %input('enter the yaxis label','s');
tt = 'Line-Profile Displacements';%input('enter the title','s');
%le = 'Single'; %input('enter the legend','s');
le2 = 'Time-Lapse 1';
le3 = 'Time-Lapse 2';
xl = xlabel(xt);
yl = ylabel(yt);
%tl = title(tt);

set(xl, 'fontweight','bold','fontsize',28,'color',fcolor);
set(yl,'fontweight','bold','fontsize',28,'color',fcolor);
%leg = legend(['\color{' fcolor '}' le2],['\color{' fcolor '}' le3],'location','west','fontcolor',fcolor);
%legend boxoff
%leg.FontSize = LegendFontSize;
%set(tl,'fontweight','bold','fontsize',title_font_size)
axis([-7 0 -10 80])
xticks([-6 -4 -2 0])

title = ['\PatternVNonPatterned ' fcolor ' on ' bcolor];
savefile = [ListPath title];
export_fig(PaVsNp,savefile,'-native');