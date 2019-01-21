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

singlesT =load('Cluster TractionStats.mat','TractionStats');
singlesT2 = load('5%wt Gels TractionStats.mat','TractionStats');

%% Correlation Coefficient
sTot = singlesT.TractionStats(:,10);
nTot = singlesT.TractionStats(:,9);
aTot = singles.vqT(:,5);

[rhoShear,pS] = corr(aTot,sTot);
[rhoNormal,pN] = corr(aTot,nTot);
[rhoSN,pSN] = corr(nTot,sTot);

rhoShear2 = rhoShear^2;
rhoNormal2 = rhoNormal^2;
rhoSN2 = rhoSN^2;

lmAS = fitlm(aTot,sTot,'Intercept',false)
lmAN = fitlm(aTot,nTot,'Intercept',false)
lmNS = fitlm(nTot,sTot,'Intercept',false)




%%
for i = 1:size(colOptions,2)
    
   
    %%SHEAR
    fcolor = colOptions{1,i};
    bcolor = colOptions{2,i};
    set(0,'defaultfigurecolor',bcolor)
    
    sheararea = figure;
    
    %Plot Data
    scatter(singles.vqT(:,5),singlesT.TractionStats(:,10),50,'square',bcolor,'markerfacecolor',fcolor)
    hold on
    
    plot([0:ceil(max(aTot)/1000)*1000],(lmAS.Coefficients{1,1})*[0:ceil(max(aTot)/1000)*1000],'LineStyle','--','Color',fcolor,'HandleVisibility','off','LineWidth',2)
%     text(30,.967,['R^2: ' sprintf('%.2f',lmAS.Rsquared.Ordinary(1,1))],'color',fcolor,'fontsize',18)
%     text(30,.867,['p: ' sprintf('%.2e',lmAS.Coefficients{2,4})],'color',fcolor,'fontsize',18)
    
    
    set(gca,'Color',bcolor,'LineWidth',2)
    %Axes, Text, Legends
    xlim([0 ceil(max(aTot)/1000)*1000])
    ylim([0 ceil(max(sTot)/750)*750])
    set(gca,'fontsize',AxisFontSize,'XColor',fcolor,'YColor',fcolor,'YMinorTick','on')
    ytickformat('%.1f')
    xt = 'Cell Area (\mum^{2})';% input('enter the xaxis label','s');
    yt = '\Sigma |Shear| (pN)'; %input('enter the yaxis label','s');
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
    title = ['\ShearTractionVsArea Cluster ' fcolor ' on ' bcolor];
    savefile = [ListPath title];
    export_fig(sheararea,savefile,'-native');
    %%
    normalarea = figure;
    
     %Plot Data
    scatter(singles.vqT(:,5),singlesT.TractionStats(:,9),50,'square',bcolor,'markerfacecolor',fcolor)
    hold on
    
    plot([0:ceil(max(aTot)/1000)*1000],lmAN.Coefficients{1,1}*[0:ceil(max(aTot)/1000)*1000],'LineStyle','--','Color',fcolor,'LineWidth',2)
%     text(30,.967,['R^2: ' sprintf('%.2f',lmAN.Rsquared.Ordinary(1,1))],'color',fcolor,'fontsize',18,'HorizontalAlignment','left')
%     text(30,.867,['p: ' sprintf('%.2e',lmAN.Coefficients{2,4})],'color',fcolor,'fontsize',18,'HorizontalAlignment','left')
%     
    xlim([0 ceil(max(aTot)/500)*500])
    ylim([0 ceil(max(nTot)/750)*750])
    set(gca,'Color',bcolor,'LineWidth',2)
    %Axes, Text, Legends
    set(gca,'fontsize',AxisFontSize,'XColor',fcolor,'YColor',fcolor,'YMinorTick','on')
    ytickformat('%.1f')
    xt = 'Cell Area (\mum^{2})';% input('enter the xaxis label','s');
    yt = '\Sigma |Normal| (pN)'; %input('enter the yaxis label','s');
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
    title = ['\NormalTractionVsArea Cluster ' fcolor ' on ' bcolor];
    savefile = [ListPath title];
    export_fig(normalarea,savefile,'-native');
    %%
    shearnormal = figure;
    
     
    %Plot Data
    scatter(singlesT.TractionStats(:,9)/max(singlesT.TractionStats(:,9)),singlesT.TractionStats(:,10)/max(singlesT.TractionStats(:,9)),50,'square',bcolor,'markerfacecolor',fcolor)
    hold on
    
    plot([0:1],lmNS.Coefficients{1,1}*[0:1]+(lmNS.Coefficients{1,1}/max(singlesT.TractionStats(:,9))),'LineStyle','--','Color',fcolor,'LineWidth',2)
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
    title = ['\ShearVsNormalTractions Cluster ' fcolor ' on ' bcolor];
    savefile = [ListPath title];
    export_fig(shearnormal,savefile,'-native');
    
    %%
    shearnormaltrue = figure;
    
     
    %Plot Data
    scatter(singlesT.TractionStats(:,9),singlesT.TractionStats(:,10),50,'square',bcolor,'markerfacecolor',fcolor)
    hold on

    plot([0:5000],lmNS.Coefficients{1,1}*[0:5000],'LineStyle','--','Color',fcolor,'LineWidth',2)
    
    xlim([0 ceil(max(nTot)/500)*500])
    ylim([0 ceil(max(sTot)/500)*500])
    
    set(gca,'Color',bcolor,'LineWidth',2)
    %Axes, Text, Legends
    set(gca,'fontsize',AxisFontSize,'XColor',fcolor,'YColor',fcolor,'YMinorTick','on')
    ytickformat('%.1f')
    xt = '\Sigma |Normal| (pN)';% input('enter the xaxis label','s');
    yt = '\Sigma |Shear| (pN)'; %input('enter the yaxis label','s');
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
    title = ['\ShearVsNormalTrueTractions Cluster ' fcolor ' on ' bcolor];
    savefile = [ListPath title];
    export_fig(shearnormaltrue,savefile,'-native');
    
end

