%%Bar graph of mean distances 
% Data is stored in dropbox publications folder
set(0,'defaultfigurecolor',[1 1 1])
fcolor = 'black';
dist = figure;
AxisFontSize = 28;
AxisTitleFontSize = 28;
LegendFontSize = 20;
bar(meanDisplacements(1,1:3),'FaceColor',[.35 .35 .35])
hold on
errorbar(meanDisplacements(1,1:3),meanDisplacements(2,1:3),'.','color',[0 0 0],'MarkerSize',1)

xt = 'Dimension';% input('enter the xaxis label','s');
yt = {'Center-to-Center'; 'Distance (\mum)'}; %input('enter the yaxis label','s');
tt = 'Line-Profile Displacements';%input('enter the title','s');
le = 'Shear'; %input('enter the legend','s');
le2 = 'Normal';
le3 = 'Border';
xl = xlabel(xt);
yl = ylabel(yt); 
%tl = title(tt);
set(gca,'xticklabel', {'X' 'Y' 'Z'})
set(xl, 'fontweight','bold','fontsize',AxisTitleFontSize); 
set(yl,'fontweight','bold','fontsize',AxisTitleFontSize);
set(gca,'fontsize',AxisFontSize,'XColor',fcolor,'YColor',fcolor,'ZColor',fcolor,'YMinorTick','on')

%leg = legend([p1 p2 p3],le,le2,le3,'location','northwest');
%leg.FontSize = LegendFontSize;
%set(tl,'fontweight','bold','fontsize',title_font_size)
ylim([0 6])
box off
filePath=cd;
title = '\Center-to-Center Spacing';
savefile = [filePath title];
%%
export_fig(dist,savefile,'-native');



%% Dark background version

%%Bar graph of mean distances 
% Data is stored in dropbox publications folder
set(0,'defaultfigurecolor',[0 0 0])
fcolor = 'white';
bcolor = 'black';

dist = figure;
AxisFontSize = 28;
AxisTitleFontSize = 28;
LegendFontSize = 20;
bar(meanDisplacements(1,1:3),'FaceColor',[.35 .35 .35])
hold on
errorbar(meanDisplacements(1,1:3),meanDisplacements(2,1:3),'.','color',fcolor,'MarkerSize',1)
set(gca,'Color',bcolor)
xt = 'Dimension';% input('enter the xaxis label','s');
yt = {'Center-to-Center'; 'Distance (\mum)'}; %input('enter the yaxis label','s');
tt = 'Line-Profile Displacements';%input('enter the title','s');
le = 'Shear'; %input('enter the legend','s');
le2 = 'Normal';
le3 = 'Border';
xl = xlabel(xt);
yl = ylabel(yt); 
%tl = title(tt);
set(gca,'xticklabel', {'X' 'Y' 'Z'})
set(xl, 'fontweight','bold','fontsize',AxisTitleFontSize,'color',fcolor); 
set(yl,'fontweight','bold','fontsize',AxisTitleFontSize,'color',fcolor);
set(gca,'fontsize',AxisFontSize,'XColor',fcolor,'YColor',fcolor,'ZColor',fcolor,'YMinorTick','on')

%leg = legend([p1 p2 p3],le,le2,le3,'location','northwest');
%leg.FontSize = LegendFontSize;
%set(tl,'fontweight','bold','fontsize',title_font_size)
ylim([0 6])
box off
filePath=cd;
title = '\Center-to-Center Spacing Dark';
savefile = [filePath title];
%%
export_fig(dist,savefile,'-native');