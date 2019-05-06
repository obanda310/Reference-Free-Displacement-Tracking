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

singles =load('Cluster ShearNormalStats.mat','vqT');
singles2 = load('5%wt Gels ShearNormalStats.mat','vqT');
singles.vqT([4 6 8 11],:) = [];
singles2.vqT([29 37 42 45 50],:) = [];
%% Correlation Coefficient
sTot = singles.vqT(:,1);
nTot = singles.vqT(:,2);
aTot = singles.vqT(:,5);

[rhoShear,pS] = corr(aTot,sTot/max(singles2.vqT(:,1)));
[rhoNormal,pN] = corr(aTot,nTot/max(singles2.vqT(:,2)));
[rhoSN,pSN] = corr(nTot,sTot);

rhoShear2 = rhoShear^2;
rhoNormal2 = rhoNormal^2;
rhoSN2 = rhoSN^2;

lmAS = fitlm(aTot,sTot/max(singles2.vqT(:,1)),'Intercept',false)
lmAN = fitlm(aTot,nTot/max(singles2.vqT(:,2)),'Intercept',false)
lmNS = fitlm(nTot,sTot,'Intercept',false)
anova(lmNS,'summary')





%%
for i = 1:size(colOptions,2)
    
    %%
    fcolor = colOptions{1,i};
    bcolor = colOptions{2,i};
    set(0,'defaultfigurecolor',bcolor)
    
    sheararea = figure;
    
    %Plot Data
    scatter(singles.vqT(:,5),singles.vqT(:,1)/max(singles2.vqT(:,1)),50,'square',bcolor,'markerfacecolor',fcolor)
    hold on
    
    plot([0:ceil(max(aTot)/1000)*1000],lmAS.Coefficients{1,1}*[0:ceil(max(aTot)/1000)*1000],'LineStyle','--','Color',fcolor,'HandleVisibility','off')
%     text(30,1.95,['R^2: ' sprintf('%.2f',lmAS.Rsquared.Ordinary(1,1))],'color',fcolor,'fontsize',18)
%     text(30,1.75,['p: ' sprintf('%.3f',lmAS.Coefficients{2,4})],'color',fcolor,'fontsize',18)
%     
    xlim([0 5000])
    ylim([0 2])
    set(gca,'Color',bcolor)
    %Axes, Text, Legends
    set(gca,'fontsize',AxisFontSize,'XColor',fcolor,'YColor',fcolor,'YMinorTick','on')
    ytickformat('%.1f')
    xt = 'Cluster Area (\mum^{2})';% input('enter the xaxis label','s');
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
    title = ['\ShearVsArea Cluster ' fcolor ' on ' bcolor];
    savefile = [ListPath title];
    export_fig(sheararea,savefile,'-native');
    %%
    normalarea = figure;        
    %Plot Data
    scatter(singles.vqT(:,5),singles.vqT(:,2)/max(singles2.vqT(:,2)),50,'square',bcolor,'markerfacecolor',fcolor)
    hold on
    plot([0:ceil(max(aTot)/1000)*1000],lmAN.Coefficients{1,1}*[0:ceil(max(aTot)/1000)*1000],'LineStyle','--','Color',fcolor,'HandleVisibility','off','LineWidth',2)
%     text(30,1.95,['R^2: ' sprintf('%.2f',lmAN.Rsquared.Ordinary(1,1))],'color',fcolor,'fontsize',18,'HorizontalAlignment','left')
%     text(30,1.75,['p: ' sprintf('%.3f',lmAN.Coefficients{2,4})],'color',fcolor,'fontsize',18,'HorizontalAlignment','left')


    xlim([0 5000])
    ylim([0 2])
    set(gca,'Color',bcolor)
    %Axes, Text, Legends
    set(gca,'fontsize',AxisFontSize,'XColor',fcolor,'YColor',fcolor,'YMinorTick','on')
    ytickformat('%.1f')
    xt = 'Cluster Area (\mum^{2})';% input('enter the xaxis label','s');
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
    title = ['\NormalVsArea Cluster ' fcolor ' on ' bcolor];
    savefile = [ListPath title];
    export_fig(normalarea,savefile,'-native');
    %%
    shearnormal = figure;
    
    
    %Plot Data
    scatter(singles.vqT(:,2)/max(singles2.vqT(:,2)),singles.vqT(:,1)/max(singles2.vqT(:,2)),50,'square',bcolor,'markerfacecolor',fcolor)
    hold on

    
    plot([0:2],lmNS.Coefficients{1,1}*[0:2],'LineStyle','--','Color',fcolor)
%     text(.01,3.9,['R^2: ' sprintf('%.2f',lmNS.Rsquared.Ordinary(1,1))],'color',fcolor,'fontsize',18)
%     text(.01,3.5,['p: ' sprintf('%.2e',lmNS.Coefficients{2,4})],'color',fcolor,'fontsize',18)
%     text(.01,3.1,['X-coef.: ' sprintf('%.2f',lmNS.Coefficients{2,1})],'color',fcolor,'fontsize',18)
    ylim([0 4])
    set(gca,'Color',bcolor)
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
    title = ['\ShearVsNormal Cluster ' fcolor ' on ' bcolor];
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
    set(gca,'Color',bcolor)
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
    title = ['\ShearVsNormalTrue Cluster ' fcolor ' on ' bcolor];
    savefile = [ListPath title];
    export_fig(shearnormaltrue,savefile,'-native');
    
end
