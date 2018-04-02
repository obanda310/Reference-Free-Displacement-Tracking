load('ShearNormalStatsBleb.mat');
vqB = vqT;

load('ShearNormalStatsSpread.mat');
vqS = vqT;

% clear all
close all
set(0,'defaultfigurecolor',[1 1 1])
ListPath = cd;

%CHANGE FONT SIZES HERE
AxisFontSize = 24;
AxisTitleFontSize = 24;
LegendFontSize = 20;

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

    scatter(0:15:90,vqS(:,1)/max(vqS(:,1)),50,'square',bcolor,'markerfacecolor','red')
    hold on
    scatter(0:15:90,vqB(:,1)/max(vqS(:,1)),50,'square',bcolor,'markerfacecolor','green')
    
    
    
    set(gca,'Color',bcolor)
    %Axes, Text, Legends
    set(gca,'fontsize',AxisFontSize,'XColor',fcolor,'YColor',fcolor,'YMinorTick','on')
    ytickformat('%.1f')
    xt = 'Time (min)';% input('enter the xaxis label','s');
    yt = '\Sigma |Shear| (AU)'; %input('enter the yaxis label','s');
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
    axis([0 90 0 1])
    xticks([0 15 30 45 60 75 90])
    
    %Export Image
    title = ['\ShearVsTime ' fcolor ' on ' bcolor];
    savefile = [ListPath title];
    export_fig(sheararea,savefile,'-native');
    
    
        %%
    fcolor = colOptions{1,i};
    bcolor = colOptions{2,i};
    set(0,'defaultfigurecolor',bcolor)
    
    normaltime = figure;
    
    %Plot Data
   
    scatter(0:15:90,vqS(:,2)/max(vqS(:,2)),50,'square',bcolor,'markerfacecolor','red')
    hold on
    scatter(0:15:90,vqB(:,2)/max(vqS(:,2)),50,'square',bcolor,'markerfacecolor','green')
    
    
    
    set(gca,'Color',bcolor)
    %Axes, Text, Legends
    set(gca,'fontsize',AxisFontSize,'XColor',fcolor,'YColor',fcolor,'YMinorTick','on')
    ytickformat('%.1f')
    xt = 'Time (min)';% input('enter the xaxis label','s');
    yt = '\Sigma |Normal| (AU)'; %input('enter the yaxis label','s');
    tt = 'Line-Profile Displacements';%input('enter the title','s');
    %le = 'Single'; %input('enter the legend','s');
    le2 = 'Time-Lapse 1';
    le3 = 'Time-Lapse 2';
    xl = xlabel(xt);
    yl = ylabel(yt);
    %tl = title(tt);
    
    set(xl, 'fontweight','bold','fontsize',28,'color',fcolor);
    set(yl,'fontweight','bold','fontsize',28,'color',fcolor);
    %leg = legend(['\color{' fcolor '}' le2],['\color{' fcolor '}' le3],'location','northeast','fontcolor',fcolor);
    %legend boxoff
    %leg.FontSize = LegendFontSize;
    %set(tl,'fontweight','bold','fontsize',title_font_size)
    
    axis([0 90 0 1])
    xticks([0 15 30 45 60 75 90])
    %Export Image
    title = ['\NormalVsTime ' fcolor ' on ' bcolor];
    savefile = [ListPath title];
    export_fig(normaltime,savefile,'-native');
    
            %%
    fcolor = colOptions{1,i};
    bcolor = colOptions{2,i};
    set(0,'defaultfigurecolor',bcolor)
    
    shearnormal = figure;
    
    %Plot Data

    scatter(vqS(:,2)/max(vqS(:,2)),vqS(:,1)/max(vqS(:,1)),50,'square',bcolor,'markerfacecolor','red')
    hold on
    scatter(vqB(:,2)/max(vqS(:,2)),vqB(:,1)/max(vqS(:,1)),50,'square',bcolor,'markerfacecolor','green')
    
    
    set(gca,'Color',bcolor)
    %Axes, Text, Legends
    set(gca,'fontsize',AxisFontSize,'XColor',fcolor,'YColor',fcolor,'YMinorTick','on')
    ytickformat('%.1f')
    yt = '\Sigma |Shear| (AU)';% input('enter the xaxis label','s');
    xt = '\Sigma |Normal| (AU)'; %input('enter the yaxis label','s');
    tt = 'Line-Profile Displacements';%input('enter the title','s');
    %le = 'Single'; %input('enter the legend','s');
    le2 = 'Time-Lapse 1';
    le3 = 'Time-Lapse 2';
    xl = xlabel(xt);
    yl = ylabel(yt);
    %tl = title(tt);
    
    set(xl, 'fontweight','bold','fontsize',28,'color',fcolor);
    set(yl,'fontweight','bold','fontsize',28,'color',fcolor);
    leg = legend(['\color{' fcolor '}' le2],['\color{' fcolor '}' le3],'location','southeast','fontcolor',fcolor);
    legend boxoff
    leg.FontSize = LegendFontSize;
    %set(tl,'fontweight','bold','fontsize',title_font_size)
    
    
    %Export Image
    title = ['\ShearVsNormal ' fcolor ' on ' bcolor];
    savefile = [ListPath title];
    export_fig(shearnormal,savefile,'-native');
end