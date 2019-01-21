clear all
close all
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

green = [.2 .7 .2];
singles =load('4hr Single Cell ShearNormalStats.mat','vqT');
%singles =load('ShearNormalStatsCluster.mat','vqT');
 
bleb =load('8hr Single Cells ShearNormalStats.mat','vqT');
spread =load('24hr Single Cells ShearNormalStats.mat','vqT');

%% Correlation Coefficient
sTot = cat(1,singles.vqT(:,1),bleb.vqT(:,1),spread.vqT(:,1));
nTot = cat(1,singles.vqT(:,2),bleb.vqT(:,2),spread.vqT(:,2));
aTot = cat(1,singles.vqT(:,5),bleb.vqT(:,5),spread.vqT(:,5));

[rhoShear,pS] = corr(aTot,sTot/max(singles.vqT(:,1)));
[rhoNormal,pN] = corr(aTot,nTot/max(singles.vqT(:,2)));
[rhoSN,pSN] = corr(nTot,sTot);

lmAS = fitlm(aTot,sTot/max(singles.vqT(:,1)),'linear')
lmAN = fitlm(aTot,nTot/max(singles.vqT(:,2)),'linear')
lmNS = fitlm(nTot,sTot,'linear')
anova(lmNS,'summary')


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
    try
    scatter(spread.vqT(:,5),spread.vqT(:,1)/max(singles.vqT(:,1)),50,'square',bcolor,'markerfacecolor','red')
    scatter(bleb.vqT(:,5),bleb.vqT(:,1)/max(singles.vqT(:,1)),50,'square',bcolor,'markerfacecolor',green)
    scatter(clusters.vqT(:,5),clusters.vqT(:,1)/max(clusters.vqT(:,1)),50,'square',bcolor,'markerfacecolor','y')
    catch
    end
    
    plot([0:ceil(max(aTot)/1000)*1000],lmAS.Coefficients{2,1}*[0:ceil(max(aTot)/1000)*1000]+(lmAS.Coefficients{1,1}/max(singles.vqT(:,1))),'LineStyle','--','Color',fcolor,'HandleVisibility','off','LineWidth',2)
%     text(30,.967,['R^2: ' sprintf('%.2f',lmAS.Rsquared.Ordinary(1,1))],'color',fcolor,'fontsize',18)
%     text(30,.867,['p: ' sprintf('%.2e',lmAS.Coefficients{2,4})],'color',fcolor,'fontsize',18)
    
    
    set(gca,'Color',bcolor,'LineWidth',2)
    %Axes, Text, Legends
    ylim([0 1])
    set(gca,'fontsize',AxisFontSize,'XColor',fcolor,'YColor',fcolor,'YMinorTick','on')
    ytickformat('%.1f')
    xt = 'Cell Area (\mum^{2})';% input('enter the xaxis label','s');
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
    try
    scatter(spread.vqT(:,5),spread.vqT(:,2)/max(singles.vqT(:,2)),50,'square',bcolor,'markerfacecolor','red')
    scatter(bleb.vqT(:,5),bleb.vqT(:,2)/max(singles.vqT(:,2)),50,'square',bcolor,'markerfacecolor',green)
    scatter(clusters.vqT(:,5),clusters.vqT(:,2)/max(clusters.vqT(:,2)),50,'square',bcolor,'markerfacecolor','y')
    catch
    end
    
    plot([0:ceil(max(aTot)/1000)*1000],lmAN.Coefficients{2,1}*[0:ceil(max(aTot)/1000)*1000]+(lmAN.Coefficients{1,1}/max(singles.vqT(:,2))),'LineStyle','--','Color',fcolor,'LineWidth',2)
%     text(30,.967,['R^2: ' sprintf('%.2f',lmAN.Rsquared.Ordinary(1,1))],'color',fcolor,'fontsize',18,'HorizontalAlignment','left')
%     text(30,.867,['p: ' sprintf('%.2e',lmAN.Coefficients{2,4})],'color',fcolor,'fontsize',18,'HorizontalAlignment','left')
%     
    ylim([0 1])
    set(gca,'Color',bcolor,'LineWidth',2)
    %Axes, Text, Legends
    set(gca,'fontsize',AxisFontSize,'XColor',fcolor,'YColor',fcolor,'YMinorTick','on')
    ytickformat('%.1f')
    xt = 'Cell Area (\mum^{2})';% input('enter the xaxis label','s');
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
%     leg = legend(['\color{' fcolor '}' le],['\color{' fcolor '}' le2],['\color{' fcolor '}' le3],'location','southeast','fontcolor',fcolor);
%     legend boxoff
%     leg.FontSize = LegendFontSize;
    %set(tl,'fontweight','bold','fontsize',title_font_size)
    
    
    %Export Image
    title = ['\NormalVsArea ' fcolor ' on ' bcolor];
    savefile = [ListPath title];
    export_fig(normalarea,savefile,'-native');
    %%
    shearnormal = figure;
    
    
    %Plot Data
    scatter(singles.vqT(:,2)/max(singles.vqT(:,2)),singles.vqT(:,1)/max(singles.vqT(:,2)),50,'square',bcolor,'markerfacecolor',fcolor)
    hold on
    try
    scatter(spread.vqT(:,2)/max(singles.vqT(:,2)),spread.vqT(:,1)/max(singles.vqT(:,2)),50,'square',bcolor,'markerfacecolor','red')
    scatter(bleb.vqT(:,2)/max(singles.vqT(:,2)),bleb.vqT(:,1)/max(singles.vqT(:,2)),50,'square',bcolor,'markerfacecolor',green)
    scatter(clusters.vqT(:,2)/max(singles.vqT(:,2)),clusters.vqT(:,1)/max(singles.vqT(:,2)),50,'square',bcolor,'markerfacecolor','y')
    catch
    end
    
    plot([0:1],lmNS.Coefficients{2,1}*[0:1]+(lmNS.Coefficients{1,1}/max(singles.vqT(:,2))),'LineStyle','--','Color',fcolor,'LineWidth',2)
%     text(.01,2.9,['R^2: ' sprintf('%.2f',lmNS.Rsquared.Ordinary(1,1))],'color',fcolor,'fontsize',18)
%     text(.01,2.6,['p: ' sprintf('%.2e',lmNS.Coefficients{2,4})],'color',fcolor,'fontsize',18)
%     text(.01,2.3,['X-coef.: ' sprintf('%.2f',lmNS.Coefficients{2,1})],'color',fcolor,'fontsize',18)
%     
    set(gca,'Color',bcolor,'LineWidth',2)
    %Axes, Text, Legends
    set(gca,'fontsize',AxisFontSize,'XColor',fcolor,'YColor',fcolor,'YMinorTick','on')
    ytickformat('%.1f')
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
%     leg = legend(['\color{' fcolor '}' le],['\color{' fcolor '}' le2],['\color{' fcolor '}' le3],'location','southeast','fontcolor',fcolor);
%     legend boxoff
%     leg.FontSize = LegendFontSize;
    %set(tl,'fontweight','bold','fontsize',title_font_size)
    
    
    %Export Image
    title = ['\ShearVsNormal ' fcolor ' on ' bcolor];
    savefile = [ListPath title];
    export_fig(shearnormal,savefile,'-native');
    
    %%
        shearnormaltrue = figure;
    
    
    %Plot Data
    scatter(singles.vqT(:,2),singles.vqT(:,1),50,'square',bcolor,'markerfacecolor',fcolor)
    hold on
    try
    scatter(spread.vqT(:,2),spread.vqT(:,1),50,'square',bcolor,'markerfacecolor','red')
    scatter(bleb.vqT(:,2),bleb.vqT(:,1),50,'square',bcolor,'markerfacecolor',green)
    scatter(clusters.vqT(:,2),clusters.vqT(:,1),50,'square',bcolor,'markerfacecolor','y')
    catch
    end
    set(gca,'Color',bcolor,'LineWidth',2)
    %Axes, Text, Legends
    set(gca,'fontsize',AxisFontSize,'XColor',fcolor,'YColor',fcolor,'YMinorTick','on')
    ytickformat('%.1f')
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
%     leg = legend(['\color{' fcolor '}' le],['\color{' fcolor '}' le2],['\color{' fcolor '}' le3],'location','southeast','fontcolor',fcolor);
%     legend boxoff
%     leg.FontSize = LegendFontSize;
    %set(tl,'fontweight','bold','fontsize',title_font_size)
    
    
    %Export Image
    title = ['\ShearVsNormalTrue ' fcolor ' on ' bcolor];
    savefile = [ListPath title];
    export_fig(shearnormaltrue,savefile,'-native');
    
end

t4m = mean(singles.vqT(:,1));
t4sd = std(singles.vqT(:,1));


t8m = mean(bleb.vqT(:,1));
t8sd = std(bleb.vqT(:,1));


t24m = mean(spread.vqT(:,1));
t24sd = std(spread.vqT(:,1));
%%
bxPlotData(:,1) = (singles.vqT(:,1));
bxPlotData(1:size(bleb.vqT,1),2) = (bleb.vqT(:,1));
bxPlotData(1:size(spread.vqT,1),3) = (spread.vqT(:,1));

bxPlotData(bxPlotData==0) = NaN;



ShearVsTime = figure
boxplot(bxPlotData)
set(gca,'fontsize',AxisFontSize,'XColor',fcolor,'YColor',fcolor,'YMinorTick','on')
    ytickformat('%.1f')
    xt = '\Sigma |Normal| (AU)';% input('enter the xaxis label','s');
    yt = '\Sigma |Shear| (AU)'; %input('enter the yaxis label','s');
    tt = 'Line-Profile Displacements';%input('enter the title','s');
    le = 'Single'; %input('enter the legend','s');
    le2 = 'Time-Lapse 1';
    le3 = 'Time-Lapse 2';
    %xl = xlabel(xt);
    yl = ylabel(yt);
    %tl = title(tt);
    xticklabels({'4hrs','8hrs','24hrs'})
    %set(xl, 'fontweight','bold','fontsize',28,'color',fcolor);
    set(yl,'fontweight','bold','fontsize',28,'color',fcolor);
    
        %Export Image
    title = ['\ShearVsTime ' fcolor ' on ' bcolor];
    savefile = [ListPath title];
    export_fig(ShearVsTime,savefile,'-native');
%%
bxPlotData2(:,1) = (singles.vqT(:,5));
bxPlotData2(1:size(bleb.vqT,1),2) = (bleb.vqT(:,5));
bxPlotData2(1:size(spread.vqT,1),3) = (spread.vqT(:,5));

bxPlotData2(bxPlotData2==0) = NaN;

AreaVsTime=figure
boxplot(bxPlotData2)
set(gca,'fontsize',AxisFontSize,'XColor',fcolor,'YColor',fcolor,'YMinorTick','on')
    ytickformat('%.1f')
    xt = '\Sigma |Normal| (AU)';% input('enter the xaxis label','s');
    yt = 'Cell Area (\mum^{2})'; %input('enter the yaxis label','s');
    tt = 'Line-Profile Displacements';%input('enter the title','s');
    le = 'Single'; %input('enter the legend','s');
    le2 = 'Time-Lapse 1';
    le3 = 'Time-Lapse 2';
    %xl = xlabel(xt);
    yl = ylabel(yt);
    %tl = title(tt);
    xticklabels({'4hrs','8hrs','24hrs'})
   % set(xl, 'fontweight','bold','fontsize',28,'color',fcolor);
    set(yl,'fontweight','bold','fontsize',28,'color',fcolor);
    
        %Export Image
    title = ['\AreaVsTime ' fcolor ' on ' bcolor];
    savefile = [ListPath title];
    export_fig(AreaVsTime,savefile,'-native');

