
clear all


load('Bleb ShearNormalStats.mat');
vqB = vqT;

load('Spread ShearNormalStats.mat');
vqS = vqT;

vqTBt = load('Bleb TractionStats.mat');
vqTB = vqTBt.TractionStats;


vqTSt = load('Spread TractionStats.mat');
vqTS = vqTSt.TractionStats;


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

green = [.2 .7 .2];

%singles =load('ShearNormalStats.mat','vqT')
uScale = 1000;
%%
for i = 1:size(colOptions,2)
    
    %%
    fcolor = colOptions{1,i};
    bcolor = colOptions{2,i};
    set(0,'defaultfigurecolor',bcolor)
    
    sheararea = figure;
    
    %Plot Data

    scatter(0:15:90,vqTS(:,10)/uScale,50,'square',bcolor,'markerfacecolor','red')
    hold on
    scatter(0:15:90,vqTB(:,10)/uScale,50,'square',bcolor,'markerfacecolor',green)
    
    
    ylim([0 .25])
    
    set(gca,'Color',bcolor)
    %Axes, Text, Legends
    set(gca,'fontsize',AxisFontSize,'XColor',fcolor,'YColor',fcolor,'YMinorTick','on','LineWidth',2)
    ytickformat('%.1f')
    xt = 'Time (min)';% input('enter the xaxis label','s');
    yt = '\Sigma |Shear| (\muN)'; %input('enter the yaxis label','s');
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
    axis([0 90 0 .25])
    xticks([0 15 30 45 60 75 90])
    
    %Export Image
    title = ['\ShearVsTimeTractions ' fcolor ' on ' bcolor];
    savefile = [ListPath title];
    export_fig(sheararea,savefile,'-native');
    
    
        %%
    fcolor = colOptions{1,i};
    bcolor = colOptions{2,i};
    set(0,'defaultfigurecolor',bcolor)
    
    normaltime = figure;
    
    %Plot Data
   
    scatter(0:15:90,vqTS(:,9)/uScale,50,'square',bcolor,'markerfacecolor','red')
    hold on
    scatter(0:15:90,vqTB(:,9)/uScale,50,'square',bcolor,'markerfacecolor',green)
    
    
    ylim([0 .25])
    
    set(gca,'Color',bcolor)
    %Axes, Text, Legends
    set(gca,'fontsize',AxisFontSize,'XColor',fcolor,'YColor',fcolor,'YMinorTick','on','LineWidth',2)
    ytickformat('%.1f')
    xt = 'Time (min)';% input('enter the xaxis label','s');
    yt = '\Sigma |Normal| (\muN)'; %input('enter the yaxis label','s');
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
    
    axis([0 90 0 .25])
    xticks([0 15 30 45 60 75 90])
    %Export Image
    title = ['\NormalVsTimeTractions ' fcolor ' on ' bcolor];
    savefile = [ListPath title];
    export_fig(normaltime,savefile,'-native');
    
            %%
    fcolor = colOptions{1,i};
    bcolor = colOptions{2,i};
    set(0,'defaultfigurecolor',bcolor)
    
    shearnormal = figure;
    
    %Plot Data

    scatter(vqTS(:,9)/uScale,vqTS(:,10)/uScale,50,'square',bcolor,'markerfacecolor','red')
    hold on
    scatter(vqTB(:,9)/uScale,vqTB(:,10)/uScale,50,'square',bcolor,'markerfacecolor',green)
    
    xlim([0 .25])
    ylim([0 .25])
    set(gca,'Color',bcolor)
    %Axes, Text, Legends
    set(gca,'fontsize',AxisFontSize,'XColor',fcolor,'YColor',fcolor,'YMinorTick','on','LineWidth',2)
    ytickformat('%.1f')
    yt = '\Sigma |Shear| (\muN)';% input('enter the xaxis label','s');
    xt = '\Sigma |Normal| (\muN)'; %input('enter the yaxis label','s');
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
    title = ['\ShearVsNormalTractionsTime ' fcolor ' on ' bcolor];
    savefile = [ListPath title];
    export_fig(shearnormal,savefile,'-native');
    
                %%
    fcolor = colOptions{1,i};
    bcolor = colOptions{2,i};
    set(0,'defaultfigurecolor',bcolor)
    
    shearnormaltime = figure;
    
    %Plot Data

    %Plot Data

    plot(0:15:90,vqTS(:,10)/vqTS(4,10),'Color','red')
    hold on
    plot(0:15:90,vqTB(:,10)/vqTB(4,10),'Color',green)
    
    plot(0:15:90,vqTS(:,9)/vqTS(4,9),'Color','red','LineStyle','--')
    hold on
    plot(0:15:90,vqTB(:,9)/vqTB(4,9),'color',green,'LineStyle','--')
    
    
    
    ylim([.5 1.5])
    set(gca,'Color',bcolor)
    %Axes, Text, Legends
    set(gca,'fontsize',AxisFontSize,'XColor',fcolor,'YColor',fcolor,'YMinorTick','on','LineWidth',2)
    ytickformat('%.1f')
    yt = '\Sigma |Shear| (\muN)';% input('enter the xaxis label','s');
    xt = '\Sigma |Normal| (\muN)'; %input('enter the yaxis label','s');
    tt = 'Line-Profile Displacements';%input('enter the title','s');
    %le = 'Single'; %input('enter the legend','s');
    le2 = 'Time-Lapse 1';
    le3 = 'Time-Lapse 2';
    xl = xlabel(xt);
    yl = ylabel(yt);
    %tl = title(tt);
    
    set(xl, 'fontweight','bold','fontsize',28,'color',fcolor);
    set(yl,'fontweight','bold','fontsize',28,'color',fcolor);
%     leg = legend(['\color{' fcolor '}' le2],['\color{' fcolor '}' le3],'location','southeast','fontcolor',fcolor);
%     legend boxoff
%     leg.FontSize = LegendFontSize;
    %set(tl,'fontweight','bold','fontsize',title_font_size)
    
    
    %Export Image
    title = ['\ShearVsNormalTractionsTime ' fcolor ' on ' bcolor];
    savefile = [ListPath title];
    export_fig(shearnormaltime,savefile,'-native');
end