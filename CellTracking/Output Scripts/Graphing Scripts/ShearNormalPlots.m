% clear all
% close all
set(0,'defaultfigurecolor',[1 1 1])
ListPath = cd;

%CHANGE FONT SIZES HERE
AxisFontSize = 24;
AxisTitleFontSize = 24;
LegendFontSize = 14;

colOptions{1,1} = 'white';
colOptions{2,1} = 'black';
colOptions{1,2} = 'black';
colOptions{2,2} = 'white';


%singles =load('ShearNormalStats.mat','vqT')

%%
for i = 1:size(colOptions,2)

%%
fcolor = colOptions{1,i};
bcolor = colOptions{2,i};
set(0,'defaultfigurecolor',bcolor)

sheararea = figure;

%Plot Data
scatter(singles.vqT(:,5),singles.vqT(:,1)/max(singles.vqT(:,1)),50,'square',bcolor,'markerfacecolor',fcolor)
hold on
scatter(spread.vqT(:,5),spread.vqT(:,1)/max(singles.vqT(:,1)),50,'square',bcolor,'markerfacecolor','red')
scatter(bleb.vqT(:,5),bleb.vqT(:,1)/max(singles.vqT(:,1)),50,'square',bcolor,'markerfacecolor','green')
set(gca,'Color',bcolor)
%Axes, Text, Legends
set(gca,'fontsize',AxisFontSize,'XColor',fcolor,'YColor',fcolor)
xt = 'Cell Spread Area (\mum^{2})';% input('enter the xaxis label','s');
yt = '\Sigma |Shear| (AU)'; %input('enter the yaxis label','s');
tt = 'Line-Profile Displacements';%input('enter the title','s');
le = 'Single'; %input('enter the legend','s');
le2 = 'Time-Lapse 1';
le3 = 'Time-Lapse 2';
xl = xlabel(xt);
yl = ylabel(yt); 
%tl = title(tt);

set(xl, 'fontweight','bold','fontsize',28,'color',fcolor); 
set(yl,'fontweight','bold','fontsize',28,'color',fcolor);
leg = legend(['\color{' fcolor '}' le],['\color{' fcolor '}' le2],['\color{' fcolor '}' le3],'location','southeast','fontcolor',fcolor);
legend boxoff
leg.FontSize = LegendFontSize;
%set(tl,'fontweight','bold','fontsize',title_font_size)


%Export Image
title = ['\ShearVsArea ' fcolor ' on ' bcolor];
savefile = [ListPath title];
export_fig(sheararea,savefile,'-native');
%%
normalarea = figure;


%Plot Data
scatter(singles.vqT(:,5),singles.vqT(:,2)/max(singles.vqT(:,2)),50,'square',bcolor,'markerfacecolor',fcolor)
hold on
scatter(spread.vqT(:,5),spread.vqT(:,2)/max(singles.vqT(:,2)),50,'square',bcolor,'markerfacecolor','red')
scatter(bleb.vqT(:,5),bleb.vqT(:,2)/max(singles.vqT(:,2)),50,'square',bcolor,'markerfacecolor','green')
set(gca,'Color',bcolor)
%Axes, Text, Legends
set(gca,'fontsize',AxisFontSize,'XColor',fcolor,'YColor',fcolor)
xt = 'Cell Spread Area (\mum^{2})';% input('enter the xaxis label','s');
yt = '\Sigma |Normal| (AU)'; %input('enter the yaxis label','s');
tt = 'Line-Profile Displacements';%input('enter the title','s');
le = 'Single'; %input('enter the legend','s');
le2 = 'Time-Lapse 1';
le3 = 'Time-Lapse 2';
xl = xlabel(xt);
yl = ylabel(yt); 
%tl = title(tt);

set(xl, 'fontweight','bold','fontsize',28,'color',fcolor); 
set(yl,'fontweight','bold','fontsize',28,'color',fcolor);
leg = legend(['\color{' fcolor '}' le],['\color{' fcolor '}' le2],['\color{' fcolor '}' le3],'location','southeast','fontcolor',fcolor);
legend boxoff
leg.FontSize = LegendFontSize;
%set(tl,'fontweight','bold','fontsize',title_font_size)


%Export Image
title = ['\NormalVsArea ' fcolor ' on ' bcolor];
savefile = [ListPath title];
export_fig(normalarea,savefile,'-native');
%%
shearnormal = figure;


%Plot Data
scatter(singles.vqT(:,2)/max(singles.vqT(:,2)),singles.vqT(:,1)/max(singles.vqT(:,1)),50,'square',bcolor,'markerfacecolor',fcolor)
hold on
scatter(spread.vqT(:,2)/max(singles.vqT(:,2)),spread.vqT(:,1)/max(singles.vqT(:,1)),50,'square',bcolor,'markerfacecolor','red')
scatter(bleb.vqT(:,2)/max(singles.vqT(:,2)),bleb.vqT(:,1)/max(singles.vqT(:,1)),50,'square',bcolor,'markerfacecolor','green')
set(gca,'Color',bcolor)
%Axes, Text, Legends
set(gca,'fontsize',AxisFontSize,'XColor',fcolor,'YColor',fcolor)
xt = '\Sigma |Normal| (AU)';% input('enter the xaxis label','s');
yt = '\Sigma |Shear| (AU)'; %input('enter the yaxis label','s');
tt = 'Line-Profile Displacements';%input('enter the title','s');
le = 'Single'; %input('enter the legend','s');
le2 = 'Time-Lapse 1';
le3 = 'Time-Lapse 2';
xl = xlabel(xt);
yl = ylabel(yt); 
%tl = title(tt);

set(xl, 'fontweight','bold','fontsize',28,'color',fcolor); 
set(yl,'fontweight','bold','fontsize',28,'color',fcolor);
leg = legend(['\color{' fcolor '}' le],['\color{' fcolor '}' le2],['\color{' fcolor '}' le3],'location','southeast','fontcolor',fcolor);
legend boxoff
leg.FontSize = LegendFontSize;
%set(tl,'fontweight','bold','fontsize',title_font_size)


%Export Image
title = ['\ShearVsNormal ' fcolor ' on ' bcolor];
savefile = [ListPath title];
export_fig(shearnormal,savefile,'-native');
end