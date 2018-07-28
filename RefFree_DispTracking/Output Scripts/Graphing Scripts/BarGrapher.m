function BarGrapher(Ys)
close all
AxisFontSize = 24;
AxisTitleFontSize = 24;
LegendFontSize = 14;
colOptions{1,1} = 'white';
colOptions{2,1} = 'black';
colOptions{1,2} = 'black';
colOptions{2,2} = 'white';

for k = 1:size(colOptions,2)
    fcolor = colOptions{1,k};
    bcolor = colOptions{2,k};
    Imbalance = figure;
    set(gcf,'unit','pixels','position',[200,200,1600,500])
    
    hold on
    
    
    bar(Ys(2,:),'FaceColor','red')
    bar(Ys(3,:),'FaceColor','blue')
    
    xt = 'Individual Cells';% input('enter the xaxis label','s');
    yt = {'\color{red}\Sigma|Displacements| (\mum)';'\color{blue}\SigmaDisplacements (\mum)'}; %input('enter the yaxis label','s');
    label{1} = xlabel(xt);
    label{2} = ylabel(yt);
    set(gca,'YMinorTick','on')
    ytickformat('%.1f')
    %errorbar(meanDisplacements(1,1:3),meanDisplacements(2,1:3),'.','color',[0 0 0],'MarkerSize',1)
    axis([0 size(Ys,2)+1 0 max(max(Ys(2,:)))])
    
    
    le{1} = 'Totals (\mum)'; %input('enter the legend','s');
    le{2} = 'Imbalance (\mum)';
    
    
    
    ColorScheme(fcolor,bcolor,label,le,AxisFontSize,LegendFontSize,1,0)
    
    %Export Image
    mkdir 'Histograms'
    title = ['\Imbalance Bar Graph ' fcolor ' on ' bcolor];
    savefile = [cd '\Histograms' title];
    export_fig(Imbalance,savefile,'-native');
end


for k = 1:size(colOptions,2)
    fcolor = colOptions{1,k};
    bcolor = colOptions{2,k};
    ImbalancePct = figure;
    set(gcf,'unit','pixels','position',[200,200,1600,500])
    
    
    
    hold on
    bar((Ys(2,:)./Ys(2,:))*100,'FaceColor','red')
    bar((Ys(3,:)./Ys(2,:))*100,'FaceColor','blue')
    %errorbar(meanDisplacements(1,1:3),meanDisplacements(2,1:3),'.','color',[0 0 0],'MarkerSize',1)
    axis([0 size(Ys,2)+1 0 100])
    xt = 'Individual Cells';% input('enter the xaxis label','s');
    yt = {'\color{red}\Sigma|Displacements| (%)';'\color{blue}\SigmaDisplacements (%)'}; %input('enter the yaxis label','s');
    
    le{1} = 'Totals (%)'; %input('enter the legend','s');
    le{2} = 'Imbalance (%)';
    
    label{1} = xlabel(xt);
    label{2} = ylabel(yt);
    set(gca,'YMinorTick','on')
    ColorScheme(fcolor,bcolor,label,le,AxisFontSize,14,1,0)
    
    %Export Image
    mkdir 'Histograms'
    title = ['\Imbalance Pct Bar Graph ' fcolor ' on ' bcolor];
    savefile = [cd '\Histograms' title];
    export_fig(ImbalancePct,savefile,'-native');
end