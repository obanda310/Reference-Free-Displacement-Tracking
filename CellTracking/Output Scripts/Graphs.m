clear all
close all

ListPath = cd;

%CHANGE FONT SIZES HERE
AxisFontSize = 24;
AxisTitleFontSize = 24;
LegendFontSize = 20;


load('ShearNormalArea.mat')

Zsum(Zsum==0)=NaN;
sheararea = figure;
%Plot Data
scatter(Area,XYsum/max(XYsum),50,'square','black','markerfacecolor','black')

%Axes, Text, Legends
set(gca,'fontsize',AxisFontSize)
xt = 'Cell Spread Area (\mum^{2})';% input('enter the xaxis label','s');
yt = '\Sigma |Shear Disp.| (AU)'; %input('enter the yaxis label','s');
tt = 'Line-Profile Displacements';%input('enter the title','s');
le = 'Shear'; %input('enter the legend','s');
le2 = 'Normal';
le3 = 'Border';
xl = xlabel(xt);
yl = ylabel(yt); 
%tl = title(tt);

set(xl, 'fontweight','bold','fontsize',AxisTitleFontSize); 
set(yl,'fontweight','bold','fontsize',AxisTitleFontSize);
%leg = legend([p1 p2 p3],le,le2,le3,'location','northwest');
%leg.FontSize = LegendFontSize;
%set(tl,'fontweight','bold','fontsize',title_font_size)


%Export Image
title = '\ShearVsArea';
savefile = [ListPath title];
export_fig(sheararea,savefile,'-native');

%%
normalarea = figure;

%Plot Data
scatter(Area,Zsum/max(Zsum),50,'square','black','markerfacecolor','black')

%Axes, Text, Legends
set(gca,'fontsize',AxisFontSize)
xt = 'Cell Spread Area (\mum^{2})';% input('enter the xaxis label','s');
yt = '\Sigma |Normal Disp.| (AU)'; %input('enter the yaxis label','s');
tt = 'Line-Profile Displacements';%input('enter the title','s');
le = 'Shear'; %input('enter the legend','s');
le2 = 'Normal';
le3 = 'Border';
xl = xlabel(xt);
yl = ylabel(yt); 
%tl = title(tt);

set(xl, 'fontweight','bold','fontsize',AxisTitleFontSize); 
set(yl,'fontweight','bold','fontsize',AxisTitleFontSize);
%leg = legend([p1 p2 p3],le,le2,le3,'location','northwest');
%leg.FontSize = LegendFontSize;
%set(tl,'fontweight','bold','fontsize',title_font_size)


%Export Image
title = '\NormalVsArea';
savefile = [ListPath title];
export_fig(normalarea,savefile,'-native');
%%
shearnormal = figure;


%Plot Data
scatter(Zsum/max(Zsum),XYsum/max(XYsum),50,'square','black','markerfacecolor','black')

%Axes, Text, Legends
set(gca,'fontsize',AxisFontSize)
xt = '\Sigma |Normal| (AU)';% input('enter the xaxis label','s');
yt = '\Sigma |Shear| (AU)'; %input('enter the yaxis label','s');
tt = 'Line-Profile Displacements';%input('enter the title','s');
le = 'Shear'; %input('enter the legend','s');
le2 = 'Normal';
le3 = 'Border';
xl = xlabel(xt);
yl = ylabel(yt); 
%tl = title(tt);

set(xl, 'fontweight','bold','fontsize',28); 
set(yl,'fontweight','bold','fontsize',28);
%leg = legend([p1 p2 p3],le,le2,le3,'location','northwest');
%leg.FontSize = LegendFontSize;
%set(tl,'fontweight','bold','fontsize',title_font_size)


%Export Image
title = '\ShearVsNormal';
savefile = [ListPath title];
export_fig(shearnormal,savefile,'-native');