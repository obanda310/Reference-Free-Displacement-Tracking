%% Methods script for creating figures to accompany discriptions of reference generation and displacement calculations
load('3Ddata.mat')
mkdir('MethodsExpFigs')
colOptions{1,1} = 'white';
colOptions{2,1} = 'black';
colOptions{1,2} = 'black';
colOptions{2,2} = 'white';
fSAxisTitles = 26;
fSAxis = 20;
filePath = cd;
sizeY = size(image.RawStack,2)*raw.dataKey(9,1);

%%
close all
% [~,figBounds] =  imcrop(HeatMap3(:,:,:,1));
% figLimits(1,1) = figBounds(1,1) *raw.dataKey(9,1);
% figLimits(1,2) = (figBounds(1,1) + figBounds(1,3))*raw.dataKey(9,1);
% figLimits(1,3) = figBounds(1,2) *raw.dataKey(9,1);
% figLimits(1,4) = (figBounds(1,2) + figBounds(1,4))*raw.dataKey(9,1);

figBounds = 0;
figLimits(1,1) = 23;
figLimits(1,2) = 26;
figLimits(1,3) = sizeY-58;
figLimits(1,4) = sizeY-55;


%save('figBounds.mat','figBounds','figLimits')

%%
%load('figBounds.mat')
figLimits(1,5) = figLimits(1,2) - figLimits(1,1);
figLimits(1,6) = figLimits(1,4) - figLimits(1,3);
figLimits(1,7) = figLimits(1,1);
figLimits(1,8) = size(image.RawStack,2)*raw.dataKey(9,1) - figLimits(1,3);
figLimits2(1,1) = sizeY - figLimits(1,4);
figLimits2(1,2) = sizeY - figLimits(1,3);

figLimits(1,3:4) = figLimits2(1,1:2);


    %% Identify Markers within figLimits
    clear regionMarks
    regionMarks = 1;
    for i = 1:size(shear.rawX,2)
        
        if min(shear.rawX(shear.rawX(:,i)>0,i))>figLimits(1,1) && min(shear.rawX(shear.rawX(:,i)>0,i))<figLimits(1,2) && min(sizeY -shear.rawY(shear.rawY(:,i)>0,i))>figLimits(1,3) && min(sizeY -shear.rawY(shear.rawY(:,i)>0,i))<figLimits(1,4)
            if max(shear.rawX(shear.rawX(:,i)>0,i))>figLimits(1,1) && max(shear.rawX(shear.rawX(:,i)>0,i))<figLimits(1,2) && max(sizeY -shear.rawY(shear.rawY(:,i)>0,i))>figLimits(1,3) && max(sizeY -shear.rawY(shear.rawY(:,i)>0,i))<figLimits(1,4)
                regionMarks = cat(1,i,regionMarks);
            end
        end
    end
    %% Create graphs with single column
    regionMarks(2:end,:) = [];
    figLimits(1,1) = min(shear.rawX(shear.rawX(:,regionMarks)>0,regionMarks))-.5;
    figLimits(1,2) = max(shear.rawX(shear.rawX(:,regionMarks)>0,regionMarks))+.5;
    figLimits(1,3) = sizeY - max(shear.rawY(shear.rawY(:,regionMarks)>0,regionMarks))-.5;
    figLimits(1,4) = sizeY - min(shear.rawY(shear.rawY(:,regionMarks)>0,regionMarks))+.5;
    

%% Create vectors to build box around cropped area
boxC = [.5 .7 .9];

    for i = [2]%:size(plane.final,2)
        bots(i,1) = min(r.r(r.r(plane.final(1:nnz(plane.final(:,i)),i),3)>0,3));
    end

    RXtb = [figLimits(1,1) figLimits(1,2) figLimits(1,2) figLimits(1,1) figLimits(1,1)];
    RXfb = [figLimits(1,1) figLimits(1,1) figLimits(1,2) figLimits(1,2) figLimits(1,1)];
    RXs1 = [figLimits(1,1) figLimits(1,1) figLimits(1,1) figLimits(1,1) figLimits(1,1)];
    RXs2 = [figLimits(1,2) figLimits(1,2) figLimits(1,2) figLimits(1,2) figLimits(1,2)];
    
    RYtb = [figLimits(1,3) figLimits(1,3) figLimits(1,4) figLimits(1,4) figLimits(1,3)];
    RYf = [figLimits(1,3) figLimits(1,3) figLimits(1,3) figLimits(1,3) figLimits(1,3)];
    RYb = [figLimits(1,4) figLimits(1,4) figLimits(1,4) figLimits(1,4) figLimits(1,4)];
    RYss = [figLimits(1,3) figLimits(1,3) figLimits(1,4) figLimits(1,4) figLimits(1,3)];
    
    
    
    bot = min(bots)
    top = max(r.r(:,3))
    RZsb = [bot bot bot bot bot];
    RZst = [top top top top top];
    RZss = [bot top top bot bot];

%%
psizeX = 400;
psizeY = 800;
scatterSize = 400;
scatterSize3D = 400;
%%
for c =1
    fcolor = colOptions{1,c};
    bcolor = colOptions{2,c};
    set(0,'defaultfigurecolor',bcolor)

    
    %%
    f0 = figure;
    set(gcf,'position',[50,50,psizeX,psizeY])
    hold on
    plot3(RXtb,RYtb,[0 0 0 0 0],'LineStyle','--','Color',boxC,'LineWidth',4)
    for i = 1:size(regionMarks,1)
        %plot3(ones(shear.numFrames,1)*shear.rawX1(i),size(res,2)*raw.dataKey(9,1)-ones(shear.numFrames,1)*shear.rawY1(i),((1:1:shear.numFrames)*raw.dataKey(10,1))','r')
        scatter3(shear.rawX(:,regionMarks(i,1)),size(res,2)*raw.dataKey(9,1)-shear.rawY(:,regionMarks(i,1)),shear.rawZ(:,regionMarks(i,1)),scatterSize,'.','MarkerEdgeColor',[0.9290    0.6940    0.1250])
    end
    
    xlim([figLimits(1,1) figLimits(1,2)])
    ylim([figLimits(1,3) figLimits(1,4)])
    xl = xlabel('X');
    yl = ylabel('Y');
    zl = zlabel('Z (\mum)');
    az = -15;
    ele = 15;
    view([az ele])
    set(gca,'Color',bcolor)
    set(gca,'XColor',fcolor,'YColor',fcolor,'ZColor',fcolor,'YMinorTick','off')
    set(xl, 'fontweight','bold','color',fcolor);
    set(yl,'fontweight','bold','color',fcolor); set(zl,'fontweight','bold','color',fcolor);
    set(gca,'fontsize',fSAxis,'LineWidth',2,'XMinorTick','off')
    set(xl,'fontsize',fSAxisTitles)
    set(yl,'fontsize',fSAxisTitles)
    set(zl,'fontsize',fSAxisTitles)
    %set(gca,'XtickLabel',[24 25 26])
    
    savefile = [strcat(filePath,'\MethodsExpFigs\')  fcolor ' on ' bcolor 'fig0.tif'];
    export_fig(f0,savefile,'-native');
    %%
    f0b = figure;
    set(gcf,'position',[50,50,psizeX,psizeY])
    hold on
    plot3(RXtb,RYtb,[0 0 0 0 0],'LineStyle','--','Color',boxC,'LineWidth',4)
    %plot3(ones(shear.numFrames,1)*shear.rawX1(i),size(res,2)*raw.dataKey(9,1)-ones(shear.numFrames,1)*shear.rawY1(i),((1:1:shear.numFrames)*raw.dataKey(10,1))','r')
    scatter3(shear.rawX(:),size(res,2)*raw.dataKey(9,1)-shear.rawY(:),shear.rawZ(:),scatterSize,'.','MarkerEdgeColor',[0.9290    0.6940    0.1250])
    
    xlim([figLimits(1,1) figLimits(1,2)])
    ylim([figLimits(1,3) figLimits(1,4)])
    xl = xlabel('X');
    yl = ylabel('Y');
    zl = zlabel('Z (\mum)');
    az = -15;
    ele = 15;
    view([az ele])
    set(gca,'Color',bcolor)
    set(gca,'XColor',fcolor,'YColor',fcolor,'ZColor',fcolor)
    set(xl, 'fontweight','bold','color',fcolor);
    set(yl,'fontweight','bold','color',fcolor); set(zl,'fontweight','bold','color',fcolor);
    set(gca,'fontsize',fSAxis,'LineWidth',2)
    set(xl,'fontsize',fSAxisTitles)
    set(yl,'fontsize',fSAxisTitles)
    set(zl,'fontsize',fSAxisTitles)
    
    savefile = [strcat(filePath,'\MethodsExpFigs\')  fcolor ' on ' bcolor 'fig0b.tif'];
    export_fig(f0b,savefile,'-native');
    
    %%
    
    f1 = figure;
    set(gcf,'position',[50,50,psizeX,psizeY])
    
    hold on
    plot3(RXtb,RYtb,[0 0 0 0 0],'LineStyle','--','Color',boxC,'LineWidth',4)
    for i = 1:size(regionMarks,1)
        plot3(ones(shear.numFrames,1)*shear.rawX1(regionMarks(i,1)),size(res,2)*raw.dataKey(9,1)-ones(shear.numFrames,1)*shear.rawY1(regionMarks(i,1)),((1:1:shear.numFrames)*raw.dataKey(10,1))','LineWidth',2,'Color',fcolor)
        scatter3(shear.rawX(:,regionMarks(i,1)),size(res,2)*raw.dataKey(9,1)-shear.rawY(:,regionMarks(i,1)),shear.rawZ(:,regionMarks(i,1)),scatterSize,'.','MarkerEdgeColor',[0.9290    0.6940    0.1250])
    end

    xlim([figLimits(1,1) figLimits(1,2)])
    ylim([figLimits(1,3) figLimits(1,4)])
    zlim([0 12])
    xl = xlabel('X');
    yl = ylabel('Y');
    zl = zlabel('Z (\mum)');
    az = -15;
    ele = 15;
    view([az ele])
    set(gca,'Color',bcolor)
    set(gca,'XColor',fcolor,'YColor',fcolor,'ZColor',fcolor)
    set(xl, 'fontweight','bold','color',fcolor);
    set(yl,'fontweight','bold','color',fcolor); set(zl,'fontweight','bold','color',fcolor);
    set(gca,'fontsize',fSAxis,'LineWidth',2)
    set(xl,'fontsize',fSAxisTitles)
    set(yl,'fontsize',fSAxisTitles)
    set(zl,'fontsize',fSAxisTitles)
    
    savefile = [strcat(filePath,'\MethodsExpFigs\')  fcolor ' on ' bcolor 'fig1.tif'];
    export_fig(f1,savefile,'-native');
    %%
    vX(:,:) = shear.rawX(:,:) - shear.ltdX(:,:);
    vY(:,:) = shear.rawY(:,:) - shear.ltdY(:,:);
    vX(vX == 0) = NaN;
    vY(vY == 0) = NaN;
    %Figure 2 - Tilt Corrected
    f2 = figure;
    set(gcf,'position',[50,50,psizeX,psizeY])
    hold on
    plot3(RXtb,RYtb,[0 0 0 0 0],'LineStyle','--','Color',boxC,'LineWidth',4)
    for i = 1:size(regionMarks,1)
        plot3(vX(:,regionMarks(i,1)),size(res,2)*raw.dataKey(9,1)-vY(:,regionMarks(i,1)),((1:1:shear.numFrames)*raw.dataKey(10,1))','Color',fcolor,'LineWidth',2)
        scatter3(shear.rawX(:,regionMarks(i,1)),size(res,2)*raw.dataKey(9,1)-shear.rawY(:,regionMarks(i,1)),shear.rawZ(:,regionMarks(i,1)),scatterSize,'.','MarkerEdgeColor',[0.9290    0.6940    0.1250])
    end
    
    xlim([figLimits(1,1) figLimits(1,2)])
    ylim([figLimits(1,3) figLimits(1,4)])
    xl = xlabel('X');
    yl = ylabel('Y');
    zl = zlabel('Z (\mum)');
    
    az = -15;
    ele = 15;
    view([az ele])
    set(gca,'Color',bcolor)
    set(gca,'XColor',fcolor,'YColor',fcolor,'ZColor',fcolor)
    set(xl, 'fontweight','bold','color',fcolor);
    set(yl,'fontweight','bold','color',fcolor); set(zl,'fontweight','bold','color',fcolor);
    set(gca,'fontsize',fSAxis,'LineWidth',2)
    set(xl,'fontsize',fSAxisTitles)
    set(yl,'fontsize',fSAxisTitles)
    set(zl,'fontsize',fSAxisTitles)
    
    savefile = [strcat(filePath,'\MethodsExpFigs\')  fcolor ' on ' bcolor 'fig2.tif'];
    export_fig(f2,savefile,'-native');
    %%
    %Figure 3 - plus Quiver
    
    f3 = figure;
    set(gcf,'position',[50,50,psizeX,psizeY])
    hold on
    plot3(RXtb,RYtb,[0 0 0 0 0],'LineStyle','--','Color',boxC,'LineWidth',4)
    for i = 1:size(regionMarks,1)
       plot3(vX(:,regionMarks(i,1)),size(res,2)*raw.dataKey(9,1)-vY(:,regionMarks(i,1)),((1:1:shear.numFrames)*raw.dataKey(10,1))','Color',fcolor,'LineWidth',2)
        scatter3(shear.rawX(:,regionMarks(i,1)),size(res,2)*raw.dataKey(9,1)-shear.rawY(:,regionMarks(i,1)),(1:1:shear.numFrames)*raw.dataKey(10,1),scatterSize,'.','MarkerEdgeColor',[0.9290    0.6940    0.1250])
        quiver3(vX(:,regionMarks(i,1)),size(res,2)*raw.dataKey(9,1)-vY(:,regionMarks(i,1)),((1:1:shear.numFrames)*raw.dataKey(10,1))',shear.ltdX(:,regionMarks(i,1)),-1*shear.ltdY(:,regionMarks(i,1)),zeros(shear.numFrames,1),0,'Color',[ 0.3010    0.7450    0.9330],'LineWidth',2)
    end
     
    xlim([figLimits(1,1) figLimits(1,2)])
    ylim([figLimits(1,3) figLimits(1,4)])
    zlim([0 12])
    xl = xlabel('X');
    yl = ylabel('Y');
    zl = zlabel('Z (\mum)');
    
    az = -15;
    ele = 15;
    view([az ele])
    set(gca,'Color',bcolor)
    set(gca,'XColor',fcolor,'YColor',fcolor,'ZColor',fcolor)
    set(xl, 'fontweight','bold','color',fcolor);
    set(yl,'fontweight','bold','color',fcolor); set(zl,'fontweight','bold','color',fcolor);
    set(gca,'fontsize',fSAxis,'LineWidth',2)
    set(xl,'fontsize',fSAxisTitles)
    set(yl,'fontsize',fSAxisTitles)
    set(zl,'fontsize',fSAxisTitles)
    
    
    savefile = [strcat(filePath,'\MethodsExpFigs\')  fcolor ' on ' bcolor 'fig3.tif'];
    export_fig(f3,savefile,'-native');
    %%
    %Figure 4 - Quiver plot only
    
    f4 = figure;
    set(gcf,'position',[50,50,psizeX,psizeY])
    hold on
    plot3(RXtb,RYtb,[0 0 0 0 0],'LineStyle','--','Color',boxC,'LineWidth',4)
    for i = 1:size(regionMarks,1)
        quiver3(vX(:,regionMarks(i,1)),size(res,2)*raw.dataKey(9,1)-vY(:,regionMarks(i,1)),((1:1:shear.numFrames)*raw.dataKey(10,1))',shear.ltdX(:,regionMarks(i,1)),-1*shear.ltdY(:,regionMarks(i,1)),zeros(shear.numFrames,1),0)
    end
     
    xlim([figLimits(1,1) figLimits(1,2)])
    ylim([figLimits(1,3) figLimits(1,4)])
    xl = xlabel('X');
    yl = ylabel('Y');
    zl = zlabel('Z (\mum)');
    az = -15;
    ele = 15;
    view([az ele])
    set(gca,'Color',bcolor)
    set(gca,'XColor',fcolor,'YColor',fcolor,'ZColor',fcolor)
    set(xl, 'fontweight','bold','color',fcolor);
    set(yl,'fontweight','bold','color',fcolor); set(zl,'fontweight','bold','color',fcolor);
    set(gca,'fontsize',fSAxis,'LineWidth',2)
    set(xl,'fontsize',fSAxisTitles)
    set(yl,'fontsize',fSAxisTitles)
    set(zl,'fontsize',fSAxisTitles)
    
    
    savefile = [strcat(filePath,'\MethodsExpFigs\')  fcolor ' on ' bcolor 'fig4.tif'];
    export_fig(f4,savefile,'-native');
    
    %%
    %Change elevation to better view 3D detections
    % az = -110;
    % ele = 15;
    
    %% Gather detection rows in window
    rZone = unique(r.row(1,r.r(:,1)>figLimits(1,1)&r.r(:,1)<figLimits(1,2)&r.r(:,2)<size(res,2)*raw.dataKey(9,1)-figLimits(1,3)&r.r(:,2)>size(res,2)*raw.dataKey(9,1)-figLimits(1,4)))';
    %%
    rZonePts = find(r.r(:,1)>figLimits(1,1)&r.r(:,1)<figLimits(1,2)&r.r(:,2)<size(res,2)*raw.dataKey(9,1)-figLimits(1,3)&r.r(:,2)>size(res,2)*raw.dataKey(9,1)-figLimits(1,4));
    %%
    %Figure 5 - 3D Detections
    
    f5 = figure;
    set(gcf,'position',[50,50,psizeX,psizeY])
    hold on
    plot3(RXtb,RYtb,[0 0 0 0 0],'LineStyle','--','Color',boxC,'LineWidth',4)
    % for i = 1:regionMarks
    %
    %     scatter3(shear.rawX(:,i),size(res,2)*raw.dataKey(9,1)-shear.rawY(:,i),(1:1:shear.numFrames)*raw.dataKey(10,1),'.')
    % end
    scatter3(shear.rawX(:,regionMarks(1,1)),size(res,2)*raw.dataKey(9,1)-shear.rawY(:,regionMarks(1,1)),shear.rawZ(:,regionMarks(1,1)),scatterSize,'.','MarkerEdgeColor',[0.9290    0.6940    0.1250])
    for i =1:size(plane.final,2)
        
        scatter3(r.r(plane.final(1:nnz(plane.final(:,i)),i),1),size(res,2)*raw.dataKey(9,1)-r.r(plane.final(1:nnz(plane.final(:,i)),i),2),r.r(plane.final(1:nnz(plane.final(:,i)),i),3),scatterSize3D,'LineWidth',2,'MarkerEdgeColor',[1 1 1])
    end
    xlim([figLimits(1,1) figLimits(1,2)])
    ylim([figLimits(1,3) figLimits(1,4)])
    xl = xlabel('X');
    yl = ylabel('Y');
    zl = zlabel('Z (\mum)');
    az = -15;
    ele = 15;
    view([az ele])
    set(gca,'Color',bcolor)
    set(gca,'XColor',fcolor,'YColor',fcolor,'ZColor',fcolor)
    set(xl, 'fontweight','bold','color',fcolor);
    set(yl,'fontweight','bold','color',fcolor); set(zl,'fontweight','bold','color',fcolor);
    set(gca,'fontsize',fSAxis,'LineWidth',2)
    set(xl,'fontsize',fSAxisTitles)
    set(yl,'fontsize',fSAxisTitles)
    set(zl,'fontsize',fSAxisTitles)
    
    
    savefile = [strcat(filePath,'\MethodsExpFigs\')  fcolor ' on ' bcolor 'fig5.tif'];
    export_fig(f5,savefile,'-native');
    %%
    %Figure 6a - 3D Detections and Row Assignment
    
    f6a = figure;
    set(gcf,'position',[50,50,psizeX,psizeY])
    hold on
    colorMapFig7 = brewermap(size(rZone,1),'*spectral');
    clear temp
    for i = 1:size(rZone,1)
        
        scatter3(r.r(rows(rZone(i,1),1:nnz(rows(rZone(i,1),:))),1),(size(image.Area,1)*raw.dataKey(9,1))-r.r(rows(rZone(i,1),1:nnz(rows(rZone(i,1),:))),2),r.r(rows(rZone(i,1),1:nnz(rows(rZone(i,1),:))),3),scatterSize3D,'MarkerEdgeColor',[colorMapFig7(i,1:3)])
        temp(1:nnz(rows(rZone(i,1),:)),1:3,i) = sortrows(r.r(rows(rZone(i,1),1:nnz(rows(rZone(i,1),:))),1:3),1);
        plot3(temp(:,1,i),(size(image.Area,1)*raw.dataKey(9,1))-temp(:,2,i),temp(:,3,i),'-','Color',[colorMapFig7(i,1:3)])
        %plot3([m3.rowFits(i,1,1) m3.rowFits(i,1,2)],[(size(image.Area,1)*raw.dataKey(9,1))-m3.rowFits(i,2,1) (size(image.Area,1)*raw.dataKey(9,1))-m3.rowFits(i,2,2)],[m3.rowFits(i,3,1) m3.rowFits(i,3,2)],'Color',[colorMapFig7(i,1:3)])
    end
    scatter3(0,0,0)
    
    xlim([figLimits(1,1) figLimits(1,2)])
    ylim([figLimits(1,3) figLimits(1,4)])
    xl = xlabel('X (\mum)');
    yl = ylabel('Y (\mum)');
    zl = zlabel('Z (\mum)');
    az = -15;
    ele = 15;
    view([az ele])
    set(gca,'Color',bcolor)
    set(gca,'XColor',fcolor,'YColor',fcolor,'ZColor',fcolor,'YMinorTick','on')
    set(xl, 'fontweight','bold','color',fcolor);
    set(yl,'fontweight','bold','color',fcolor); set(zl,'fontweight','bold','color',fcolor);
    set(gca,'fontsize',fSAxis,'LineWidth',2,'XMinorTick','on')
    set(xl,'fontsize',fSAxisTitles)
    set(yl,'fontsize',fSAxisTitles)
    set(zl,'fontsize',fSAxisTitles)
    
    
    temp(temp==0) = NaN;
    
    savefile = [strcat(filePath,'\MethodsExpFigs\')  fcolor ' on ' bcolor 'fig6a.tif'];
    export_fig(f6a,savefile,'-native');
    %%
    %Figure 6b - 3D Detections and Row Fits
    
    az2=-15;
    ele2=15;
    
    f6b = figure;
    set(gcf,'position',[50,50,psizeX,psizeY])
    hold on
    colorMapFig7 = brewermap(size(rZone,1),'*spectral');
    rowsFilt = rows;
    for i = 1:size(rZone,1)
        rowsFilt(rZone(i,1),:) = 0;
        rowsFilt(rZone(i,1),1:nnz(intersect(rZonePts,rows(rZone(i,1),:)))) = intersect(rZonePts,rows(rZone(i,1),:));
    end
    
    for i = 1:size(rZone,1)
        scatter3(r.r(rows(rZone(i,1),1:nnz(rows(rZone(i,1),:))),1),(size(image.Area,1)*raw.dataKey(9,1))-r.r(rows(rZone(i,1),1:nnz(rows(rZone(i,1),:))),2),r.r(rows(rZone(i,1),1:nnz(rows(rZone(i,1),:))),3),scatterSize3D,'MarkerEdgeColor',[.5 .9 .5])
        plot3(temp(:,1,i),(size(image.Area,1)*raw.dataKey(9,1))-temp(:,2,i),temp(:,3,i),'-','Color',[.5 .9 .5])
        plot3([m3.rowFits(rZone(i,1),1,1) m3.rowFits(rZone(i,1),1,2)],[(size(image.Area,1)*raw.dataKey(9,1))-m3.rowFits(rZone(i,1),2,1) (size(image.Area,1)*raw.dataKey(9,1))-m3.rowFits(rZone(i,1),2,2)],[m3.rowFits(rZone(i,1),3,1) m3.rowFits(rZone(i,1),3,2)],'r')
    end
    
    %     for i = 1:size(rZone,1)
    %         scatter3(r.r(rowsFilt(rZone(i,1),1:nnz(rowsFilt(rZone(i,1),:))),1),(size(image.Area,1)*raw.dataKey(9,1))-r.r(rowsFilt(rZone(i,1),1:nnz(rowsFilt(rZone(i,1),:))),2),r.r(rowsFilt(rZone(i,1),1:nnz(rowsFilt(rZone(i,1),:))),3),'MarkerEdgeColor',[colorMapFig7(i,1:3)])
    %         temp2(:,1:3) = sortrows(r.r(rowsFilt(rZone(i,1),1:nnz(rowsFilt(rZone(i,1),:))),1:3),1);
    %         plot3(temp2(:,1),(size(image.Area,1)*raw.dataKey(9,1))-temp2(:,2),temp2(:,3),'Color',[colorMapFig7(i,1:3)])
    %     end
    
    xlim([figLimits(1,1) figLimits(1,2)])
    ylim([figLimits(1,3) figLimits(1,4)])
    xl = xlabel('X (\mum)');
    yl = ylabel('Y (\mum)');
    zl = zlabel('Z (\mum)');
    
    az = -15;
    ele = 15;
    view([az ele])
    set(gca,'Color',bcolor)
    set(gca,'XColor',fcolor,'YColor',fcolor,'ZColor',fcolor,'YMinorTick','on')
    set(xl, 'fontweight','bold','color',fcolor);
    set(yl,'fontweight','bold','color',fcolor); set(zl,'fontweight','bold','color',fcolor);
    set(gca,'fontsize',fSAxis,'LineWidth',2,'XMinorTick','on')
    set(xl,'fontsize',fSAxisTitles)
    set(yl,'fontsize',fSAxisTitles)
    set(zl,'fontsize',fSAxisTitles)
    
    
    savefile = [strcat(filePath,'\MethodsExpFigs\')  fcolor ' on ' bcolor 'fig6b.tif'];
    export_fig(f6b,savefile,'-native');
    %%
     %Figure 6C - 3D Detections and Row Fits
    
    f6c = figure;
    set(gcf,'position',[50,50,900,300])
    hold on
    colorMapFig7 = brewermap(size(rZone,1),'*spectral');
    rowsFilt = rows;
    plot3(RXtb,RYtb,[0 0 0 0 0],'LineStyle','--','Color',boxC,'LineWidth',2)
        for i = 1:size(regionMarks,1)
       %plot3(vX(:,regionMarks(i,1)),size(res,2)*raw.dataKey(9,1)-vY(:,regionMarks(i,1)),((1:1:shear.numFrames)*raw.dataKey(10,1))','Color',fcolor,'LineWidth',2)
        %scatter3(shear.rawX(:,regionMarks(i,1)),size(res,2)*raw.dataKey(9,1)-shear.rawY(:,regionMarks(i,1)),(1:1:shear.numFrames)*raw.dataKey(10,1),scatterSize,'.','MarkerEdgeColor',[0.9290    0.6940    0.1250])
        %quiver3(vX(:,regionMarks(i,1)),size(res,2)*raw.dataKey(9,1)-vY(:,regionMarks(i,1)),((1:1:shear.numFrames)*raw.dataKey(10,1))',shear.ltdX(:,regionMarks(i,1)),-1*shear.ltdY(:,regionMarks(i,1)),zeros(shear.numFrames,1),0,'Color',[ 0.3010    0.7450    0.9330],'LineWidth',2)
    end
    
    for i = 1:size(rZone,1)
        rowsFilt(rZone(i,1),:) = 0;
        rowsFilt(rZone(i,1),1:nnz(intersect(rZonePts,rows(rZone(i,1),:)))) = intersect(rZonePts,rows(rZone(i,1),:));
    end
    
    for i = 1:size(rZone,1)
        
        scatter3(r.r(rows(rZone(i,1),1:nnz(rows(rZone(i,1),:))),1),(size(image.Area,1)*raw.dataKey(9,1))-r.r(rows(rZone(i,1),1:nnz(rows(rZone(i,1),:))),2),r.r(rows(rZone(i,1),1:nnz(rows(rZone(i,1),:))),3),500,'.','MarkerEdgeColor',[1 1 1])
        %plot3(temp(:,1,i),(size(image.Area,1)*raw.dataKey(9,1))-temp(:,2,i),temp(:,3,i),'-','Color',[.5 .5 .5])
        %plot3([m3.rowFits(rZone(i,1),1,1) m3.rowFits(rZone(i,1),1,2)],[(size(image.Area,1)*raw.dataKey(9,1))-m3.rowFits(rZone(i,1),2,1) (size(image.Area,1)*raw.dataKey(9,1))-m3.rowFits(rZone(i,1),2,2)],[m3.rowFits(rZone(i,1),3,1) m3.rowFits(rZone(i,1),3,2)],'Color','r','LineWidth',2)
    end
    
    for i = 1:size(rZone,1)
        clear temp2
        scatter3(r.r(rowsFilt(rZone(i,1),1:nnz(rowsFilt(rZone(i,1),:))),1),(size(image.Area,1)*raw.dataKey(9,1))-r.r(rowsFilt(rZone(i,1),1:nnz(rowsFilt(rZone(i,1),:))),2),r.r(rowsFilt(rZone(i,1),1:nnz(rowsFilt(rZone(i,1),:))),3),scatterSize3D,'MarkerEdgeColor',[1 1 1],'LineWidth',3)
        temp2(:,1:3) = sortrows(r.r(rowsFilt(rZone(i,1),1:nnz(rowsFilt(rZone(i,1),:))),1:3),1);
        %plot3(temp2(:,1),(size(image.Area,1)*raw.dataKey(9,1))-temp2(:,2),temp2(:,3),'Color',[colorMapFig7(i,1:3)])
    end
    
    %xlim([figLimits(1,1) figLimits(1,2)])
    %ylim([figLimits(1,3) figLimits(1,4)])
    xlim([0 100])
    ylim([52 63])
    xl = xlabel('X (\mum)');
    yl = ylabel('Y (\mum)');
    zl = zlabel('Z (\mum)');
    
    az = -15;
    ele = 15;
    view([az ele])
    set(gca,'Color',bcolor)
    set(gca,'XColor',fcolor,'YColor',fcolor,'ZColor',fcolor)
    set(xl, 'fontweight','bold','color',fcolor);
    set(yl,'fontweight','bold','color',fcolor); set(zl,'fontweight','bold','color',fcolor);
    set(gca,'fontsize',fSAxis,'LineWidth',2)
    set(xl,'fontsize',fSAxisTitles)
    set(yl,'fontsize',fSAxisTitles)
    set(zl,'fontsize',fSAxisTitles)
    
    
    savefile = [strcat(filePath,'\MethodsExpFigs\')  fcolor ' on ' bcolor 'fig6c2.tif'];
    export_fig(f6c,savefile,'-native');
    %%
    %Figure 6C - 3D Detections and Row Fits
    
    f6c = figure;
    set(gcf,'position',[50,50,900,300])
    hold on
    colorMapFig7 = brewermap(size(rZone,1),'*spectral');
    rowsFilt = rows;
    plot3(RXtb,RYtb,[0 0 0 0 0],'LineStyle','--','Color',boxC,'LineWidth',2)
        for i = 1:size(regionMarks,1)
       plot3(vX(:,regionMarks(i,1)),size(res,2)*raw.dataKey(9,1)-vY(:,regionMarks(i,1)),((1:1:shear.numFrames)*raw.dataKey(10,1))','Color',fcolor,'LineWidth',2)
        scatter3(shear.rawX(:,regionMarks(i,1)),size(res,2)*raw.dataKey(9,1)-shear.rawY(:,regionMarks(i,1)),(1:1:shear.numFrames)*raw.dataKey(10,1),scatterSize,'.','MarkerEdgeColor',[0.9290    0.6940    0.1250])
        %quiver3(vX(:,regionMarks(i,1)),size(res,2)*raw.dataKey(9,1)-vY(:,regionMarks(i,1)),((1:1:shear.numFrames)*raw.dataKey(10,1))',shear.ltdX(:,regionMarks(i,1)),-1*shear.ltdY(:,regionMarks(i,1)),zeros(shear.numFrames,1),0,'Color',[ 0.3010    0.7450    0.9330],'LineWidth',2)
    end
    
    for i = 1:size(rZone,1)
        rowsFilt(rZone(i,1),:) = 0;
        rowsFilt(rZone(i,1),1:nnz(intersect(rZonePts,rows(rZone(i,1),:)))) = intersect(rZonePts,rows(rZone(i,1),:));
    end
    
    for i = 1:size(rZone,1)
        
        scatter3(r.r(rows(rZone(i,1),1:nnz(rows(rZone(i,1),:))),1),(size(image.Area,1)*raw.dataKey(9,1))-r.r(rows(rZone(i,1),1:nnz(rows(rZone(i,1),:))),2),r.r(rows(rZone(i,1),1:nnz(rows(rZone(i,1),:))),3),500,'.','MarkerEdgeColor',[1 1 1])
        plot3(temp(:,1,i),(size(image.Area,1)*raw.dataKey(9,1))-temp(:,2,i),temp(:,3,i),'-','Color',[.5 .5 .5])
        plot3([m3.rowFits(rZone(i,1),1,1) m3.rowFits(rZone(i,1),1,2)],[(size(image.Area,1)*raw.dataKey(9,1))-m3.rowFits(rZone(i,1),2,1) (size(image.Area,1)*raw.dataKey(9,1))-m3.rowFits(rZone(i,1),2,2)],[m3.rowFits(rZone(i,1),3,1) m3.rowFits(rZone(i,1),3,2)],'Color','r','LineWidth',2)
    end
    
    for i = 1:size(rZone,1)
        clear temp2
        scatter3(r.r(rowsFilt(rZone(i,1),1:nnz(rowsFilt(rZone(i,1),:))),1),(size(image.Area,1)*raw.dataKey(9,1))-r.r(rowsFilt(rZone(i,1),1:nnz(rowsFilt(rZone(i,1),:))),2),r.r(rowsFilt(rZone(i,1),1:nnz(rowsFilt(rZone(i,1),:))),3),scatterSize3D,'MarkerEdgeColor',[colorMapFig7(i,1:3)])
        temp2(:,1:3) = sortrows(r.r(rowsFilt(rZone(i,1),1:nnz(rowsFilt(rZone(i,1),:))),1:3),1);
        plot3(temp2(:,1),(size(image.Area,1)*raw.dataKey(9,1))-temp2(:,2),temp2(:,3),'Color',[colorMapFig7(i,1:3)])
    end
    
    %xlim([figLimits(1,1) figLimits(1,2)])
    %ylim([figLimits(1,3) figLimits(1,4)])
    xlim([0 100])
    ylim([52 63])
    xl = xlabel('X (\mum)');
    yl = ylabel('Y (\mum)');
    zl = zlabel('Z (\mum)');
    
    az = -15;
    ele = 15;
    view([az ele])
    set(gca,'Color',bcolor)
    set(gca,'XColor',fcolor,'YColor',fcolor,'ZColor',fcolor)
    set(xl, 'fontweight','bold','color',fcolor);
    set(yl,'fontweight','bold','color',fcolor); set(zl,'fontweight','bold','color',fcolor);
    set(gca,'fontsize',fSAxis,'LineWidth',2)
    set(xl,'fontsize',fSAxisTitles)
    set(yl,'fontsize',fSAxisTitles)
    set(zl,'fontsize',fSAxisTitles)
    
    
    savefile = [strcat(filePath,'\MethodsExpFigs\')  fcolor ' on ' bcolor 'fig6c.tif'];
    export_fig(f6c,savefile,'-native');
    %%
    %Figure 7a - Row Fits/Column Fits
    clear rMem
    for i = 1:size(rZone,1)
        
        rMem(i,1:size(find(r.row == rZone(i,1)),2)) = find(r.row == rZone(i,1));   
    end
    rMem2 = rMem(:);
    rMem2(rMem2(:,1)==0,:) = [];
    for i = 1:size(rMem,2)
    rMemCol(i,1) = r.col(rMem(1,i));
    end
   
    f7a = figure;
    set(gcf,'position',[50,50,900,300])
    hold on
    colorMapFig7 = brewermap(size(rows,1),'*spectral');
    
    for i = 1:size(rMemCol,1)
        if rMemCol(i,1)>0
            i
        plot3(vX(:,rMemCol(i,1)),size(res,2)*raw.dataKey(9,1)-vY(:,rMemCol(i,1)),((1:1:shear.numFrames)*raw.dataKey(10,1))','Color',[.5 .5 .5],'LineWidth',2)
        end
    end
    
    
    for i = 1:size(rZone,1)
        
        scatter3(r.r(rows(rZone(i,1),1:nnz(rows(rZone(i,1),:))),1),(size(image.Area,1)*raw.dataKey(9,1))-r.r(rows(rZone(i,1),1:nnz(rows(rZone(i,1),:))),2),r.r(rows(rZone(i,1),1:nnz(rows(rZone(i,1),:))),3),500,'.','MarkerEdgeColor',[1 1 1])
        plot3(temp(:,1,i),(size(image.Area,1)*raw.dataKey(9,1))-temp(:,2,i),temp(:,3,i),'-','Color',[.5 .5 .5])
        plot3([m3.rowFits(rZone(i,1),1,1) m3.rowFits(rZone(i,1),1,2)],[(size(image.Area,1)*raw.dataKey(9,1))-m3.rowFits(rZone(i,1),2,1) (size(image.Area,1)*raw.dataKey(9,1))-m3.rowFits(rZone(i,1),2,2)],[m3.rowFits(rZone(i,1),3,1) m3.rowFits(rZone(i,1),3,2)],'Color','r','LineWidth',2)
    end
    
    scatter3(m3.refSC(rMem2,1),(size(image.Area,1)*raw.dataKey(9,1))-m3.refSC(rMem2,2),m3.refSC(rMem2,3),100,'x','MarkerEdgeColor',[0 1 0],'LineWidth',2)
    
            for i = 1:size(regionMarks,1)
       plot3(vX(:,regionMarks(i,1)),size(res,2)*raw.dataKey(9,1)-vY(:,regionMarks(i,1)),((1:1:shear.numFrames)*raw.dataKey(10,1))','Color',fcolor,'LineWidth',2)
        scatter3(shear.rawX(:,regionMarks(i,1)),size(res,2)*raw.dataKey(9,1)-shear.rawY(:,regionMarks(i,1)),(1:1:shear.numFrames)*raw.dataKey(10,1),scatterSize,'.','MarkerEdgeColor',[0.9290    0.6940    0.1250])
        %quiver3(vX(:,regionMarks(i,1)),size(res,2)*raw.dataKey(9,1)-vY(:,regionMarks(i,1)),((1:1:shear.numFrames)*raw.dataKey(10,1))',shear.ltdX(:,regionMarks(i,1)),-1*shear.ltdY(:,regionMarks(i,1)),zeros(shear.numFrames,1),0,'Color',[ 0.3010    0.7450    0.9330],'LineWidth',2)
    end
    scatter3(0,0,0)
    
    xlim([0 100])
    ylim([52 63])
    xl = xlabel('X (\mum)');
    yl = ylabel('Y (\mum)');
    zl = zlabel('Z (\mum)');
    
    az = -15;
    ele = 15;
    view([az ele])
    set(gca,'Color',bcolor)
    set(gca,'XColor',fcolor,'YColor',fcolor,'ZColor',fcolor,'YMinorTick','on')
    set(xl, 'fontweight','bold','color',fcolor);
    set(yl,'fontweight','bold','color',fcolor); set(zl,'fontweight','bold','color',fcolor);
    set(gca,'fontsize',fSAxis,'LineWidth',2,'XMinorTick','on')
    set(xl,'fontsize',fSAxisTitles)
    set(yl,'fontsize',fSAxisTitles)
    set(zl,'fontsize',fSAxisTitles)
    %%
    savefile = [strcat(filePath,'\MethodsExpFigs\')  fcolor ' on ' bcolor 'fig7a.tif'];
    export_fig(f7a,savefile,'-native');
    
        %%
    %Figure 7a2 - Row Fits/Column Fits Single Column
    clear rMem
    for i = 1:size(rZone,1)
        
        rMem(i,1:size(find(r.row == rZone(i,1)),2)) = find(r.row == rZone(i,1));   
    end
    rMem2 = rMem(:);
    for i = 1:size(rMem,2)
    rMemCol(i,1) = r.col(rMem(1,i));
    end
   
    f7a = figure;
    set(gcf,'position',[50,50,psizeX,psizeY])
    hold on
    plot3(RXtb,RYtb,[0 0 0 0 0],'LineStyle','--','Color',boxC,'LineWidth',4)
    colorMapFig7 = brewermap(size(rows,1),'*spectral');
    
    for i = 1:size(rMemCol,1)
        if rMemCol(i,1)>0
            i
        %plot3(vX(:,rMemCol(i,1)),size(res,2)*raw.dataKey(9,1)-vY(:,rMemCol(i,1)),((1:1:shear.numFrames)*raw.dataKey(10,1))','Color',[.5 .5 .5],'LineWidth',2)
        end
    end
    
    for i = 1:size(rZone,1)
        
        scatter3(r.r(rows(rZone(i,1),1:nnz(rows(rZone(i,1),:))),1),(size(image.Area,1)*raw.dataKey(9,1))-r.r(rows(rZone(i,1),1:nnz(rows(rZone(i,1),:))),2),r.r(rows(rZone(i,1),1:nnz(rows(rZone(i,1),:))),3),500,'o','MarkerEdgeColor',[1 1 1],'LineWidth',2)
        %plot3(temp(:,1,i),(size(image.Area,1)*raw.dataKey(9,1))-temp(:,2,i),temp(:,3,i),'-','Color',[.5 .5 .5])
        plot3([m3.rowFits(rZone(i,1),1,1) m3.rowFits(rZone(i,1),1,2)],[(size(image.Area,1)*raw.dataKey(9,1))-m3.rowFits(rZone(i,1),2,1) (size(image.Area,1)*raw.dataKey(9,1))-m3.rowFits(rZone(i,1),2,2)],[m3.rowFits(rZone(i,1),3,1) m3.rowFits(rZone(i,1),3,2)],'Color','r','LineWidth',2)
    end
    
    scatter3(m3.refSC(rZonePts,1),(size(image.Area,1)*raw.dataKey(9,1))-m3.refSC(rZonePts,2),m3.refSC(rZonePts,3),800,'x','MarkerEdgeColor',[0 1 0],'LineWidth',3)
    
            for i = 1:size(regionMarks,1)
       plot3(vX(:,regionMarks(i,1)),size(res,2)*raw.dataKey(9,1)-vY(:,regionMarks(i,1)),((1:1:shear.numFrames)*raw.dataKey(10,1))','Color',fcolor,'LineWidth',2)
        scatter3(shear.rawX(:,regionMarks(i,1)),size(res,2)*raw.dataKey(9,1)-shear.rawY(:,regionMarks(i,1)),(1:1:shear.numFrames)*raw.dataKey(10,1),scatterSize,'.','MarkerEdgeColor',[0.9290    0.6940    0.1250])
        %quiver3(vX(:,regionMarks(i,1)),size(res,2)*raw.dataKey(9,1)-vY(:,regionMarks(i,1)),((1:1:shear.numFrames)*raw.dataKey(10,1))',shear.ltdX(:,regionMarks(i,1)),-1*shear.ltdY(:,regionMarks(i,1)),zeros(shear.numFrames,1),0,'Color',[ 0.3010    0.7450    0.9330],'LineWidth',2)
    end
    scatter3(0,0,0)
    
     xlim([figLimits(1,1) figLimits(1,2)])
    ylim([figLimits(1,3) figLimits(1,4)])
    zlim([0 12])
    xl = xlabel('X');
    yl = ylabel('Y');
    zl = zlabel('Z (\mum)');
    
    az = -15;
    ele = 15;
    view([az ele])
    set(gca,'Color',bcolor)
    set(gca,'XColor',fcolor,'YColor',fcolor,'ZColor',fcolor)
    set(xl, 'fontweight','bold','color',fcolor);
    set(yl,'fontweight','bold','color',fcolor); set(zl,'fontweight','bold','color',fcolor);
    set(gca,'fontsize',fSAxis,'LineWidth',2)
    set(xl,'fontsize',fSAxisTitles)
    set(yl,'fontsize',fSAxisTitles)
    set(zl,'fontsize',fSAxisTitles)
    
    savefile = [strcat(filePath,'\MethodsExpFigs\')  fcolor ' on ' bcolor 'fig7a2.tif'];
    export_fig(f7a,savefile,'-native');
    
    
    %%
    %Figure 7b - 3D Detections/Row Fits/Column Fits
    
    f7b = figure;
    set(gcf,'position',[50,50,psizeX,psizeY])
    hold on
    colorMapFig7 = brewermap(size(rows,1),'*spectral');
    
    for i = 1:size(m3.rowFits,1)
        scatter3(r.r(rows(i,1:nnz(rows(i,:))),1),(size(image.Area,1)*raw.dataKey(9,1))-r.r(rows(i,1:nnz(rows(i,:))),2),r.r(rows(i,1:nnz(rows(i,:))),3),scatterSize3D,'MarkerEdgeColor',[colorMapFig7(i,1:3)])
        plot3([m3.rowFits(i,1,1) m3.rowFits(i,1,2)],[(size(image.Area,1)*raw.dataKey(9,1))-m3.rowFits(i,2,1) (size(image.Area,1)*raw.dataKey(9,1))-m3.rowFits(i,2,2)],[m3.rowFits(i,3,1) m3.rowFits(i,3,2)],'Color',[colorMapFig7(i,1:3)])
    end
    scatter3(m3.refSC(:,1),(size(image.Area,1)*raw.dataKey(9,1))-m3.refSC(:,2),m3.refSC(:,3),scatterSize3D,'x','r')
    scatter3(0,0,0)
    
    xlim([figLimits(1,1) figLimits(1,2)])
    ylim([figLimits(1,3) figLimits(1,4)])
    xl = xlabel('X (\mum)');
    yl = ylabel('Y (\mum)');
    zl = zlabel('Z (\mum)');
    
    az = -15;
    ele = 15;
    view([az ele])
    set(gca,'Color',bcolor)
    set(gca,'XColor',fcolor,'YColor',fcolor,'ZColor',fcolor)
    set(xl, 'fontweight','bold','color',fcolor);
    set(yl,'fontweight','bold','color',fcolor); set(zl,'fontweight','bold','color',fcolor);
    set(gca,'fontsize',fSAxis,'LineWidth',2)
    set(xl,'fontsize',fSAxisTitles)
    set(yl,'fontsize',fSAxisTitles)
    set(zl,'fontsize',fSAxisTitles)
    
    savefile = [strcat(filePath,'\MethodsExpFigs\')   fcolor ' on ' bcolor 'fig7b.tif'];
    export_fig(f7b,savefile,'-native');
    %%
    %Figure 8 - 3D Detections/Row Fits/Quiver
    f8 = figure;
    set(gcf,'position',[50,50,psizeX,psizeY])
    hold on
    plot3(RXtb,RYtb,[0 0 0 0 0],'LineStyle','--','Color',boxC,'LineWidth',4)
    colorMapFig7 = brewermap(size(rows,1),'*spectral');
    
    for i = 1:size(m3.rowFits,1)
        scatter3(r.r(rows(i,1:nnz(rows(i,:))),1),(size(image.Area,1)*raw.dataKey(9,1))-r.r(rows(i,1:nnz(rows(i,:))),2),r.r(rows(i,1:nnz(rows(i,:))),3),400,'MarkerEdgeColor',[1 1 1],'LineWidth',2)
        %plot3([m3.rowFits(i,1,1) m3.rowFits(i,1,2)],[(size(image.Area,1)*raw.dataKey(9,1))-m3.rowFits(i,2,1) (size(image.Area,1)*raw.dataKey(9,1))-m3.rowFits(i,2,2)],[m3.rowFits(i,3,1) m3.rowFits(i,3,2)],'Color',[colorMapFig7(i,1:3)])
    end
    
    scatter3(m3.refSC(rZonePts,1),(size(image.Area,1)*raw.dataKey(9,1))-m3.refSC(rZonePts,2),m3.refSC(rZonePts,3),800,'x','g','LineWidth',3)
    quiver3(m3.refSC(rZonePts,1),(size(image.Area,1)*raw.dataKey(9,1))-m3.refSC(rZonePts,2),m3.refSC(rZonePts,3),m3.disp(rZonePts,1),m3.disp(rZonePts,2)*-1,m3.disp(rZonePts,3),0,fcolor,'LineWidth',3)
     xlim([figLimits(1,1) figLimits(1,2)])
    ylim([figLimits(1,3) figLimits(1,4)])
%     zlim([3 12])
    xl = xlabel('X');
    yl = ylabel('Y');
    zl = zlabel('Z (\mum)');
    plot3(RXtb,RYtb,[0 0 0 0 0],'LineStyle','--','Color',boxC,'LineWidth',4)
    az = -15;
    ele = 15;
    view([az ele])
    set(gca,'Color',bcolor)
    set(gca,'XColor',fcolor,'YColor',fcolor,'ZColor',fcolor)
    set(xl, 'fontweight','bold','color',fcolor);
    set(yl,'fontweight','bold','color',fcolor); set(zl,'fontweight','bold','color',fcolor);
    set(gca,'fontsize',fSAxis,'LineWidth',2)
    set(xl,'fontsize',fSAxisTitles)
    set(yl,'fontsize',fSAxisTitles)
    set(zl,'fontsize',fSAxisTitles)
    
    
    savefile = [strcat(filePath,'\MethodsExpFigs\')  fcolor ' on ' bcolor 'fig8.tif'];
    export_fig(f8,savefile,'-native');
    
    %%
    %Figure 9 - Quiver Plot Full
    az= -15;
    ele= 50;
    
    f9 = figure;
    hold on
    for i = 1:size(plane.final,2)
        quiver3(m3.refSC(plane.final(1:nnz(plane.final(:,i)),i),1),(size(image.Area,1)*raw.dataKey(9,1))-m3.refSC(plane.final(1:nnz(plane.final(:,i)),i),2),m3.refSC(plane.final(1:nnz(plane.final(:,i)),i),3),m3.disp(plane.final(1:nnz(plane.final(:,i)),i),1),m3.disp(plane.final(1:nnz(plane.final(:,i)),i),2)*-1,m3.disp(plane.final(1:nnz(plane.final(:,i)),i),3),0)
    end
    xl = xlabel('X (\mum)');
    yl = ylabel('Y (\mum)');
    zl = zlabel('Z (\mum)');
    view([az ele])
    set(gca,'Color',bcolor)
    set(gca,'XColor',fcolor,'YColor',fcolor,'ZColor',fcolor,'YMinorTick','on')
    set(xl, 'fontweight','bold','color',fcolor);
    set(yl,'fontweight','bold','color',fcolor); set(zl,'fontweight','bold','color',fcolor);
    set(gca,'fontsize',fSAxis,'LineWidth',2,'XMinorTick','on')
    set(xl,'fontsize',fSAxisTitles)
    set(yl,'fontsize',fSAxisTitles)
    set(zl,'fontsize',fSAxisTitles)
    
    savefile = [strcat(filePath,'\MethodsExpFigs\')  fcolor ' on ' bcolor 'fig9.tif'];
    export_fig(f9,savefile,'-native');
    
    %%
    %Figure 10a - Quiver Plot Full With Profile plane Trace
    
    az= -12;
    ele= 26;
    
    
    
    f9 = figure;
    
    hold on
    for i = 1:size(plane.final,2)
        quiver3(m3.refSC(plane.final(1:nnz(plane.final(:,i)),i),1),(size(image.Area,1)*raw.dataKey(9,1))-m3.refSC(plane.final(1:nnz(plane.final(:,i)),i),2),m3.refSC(plane.final(1:nnz(plane.final(:,i)),i),3),m3.disp(plane.final(1:nnz(plane.final(:,i)),i),1),m3.disp(plane.final(1:nnz(plane.final(:,i)),i),2)*-1,m3.disp(plane.final(1:nnz(plane.final(:,i)),i),3),0)
        bots(i,1) = min(r.r(r.r(plane.final(1:nnz(plane.final(:,i)),i),3)>0,3));
    end
    xl = xlabel('X (\mum)');
    yl = ylabel('Y');
    zl = zlabel('Z (\mum)');
    view([az ele])
    set(gca,'Color',bcolor)
    set(gca,'XColor',fcolor,'YColor',fcolor,'ZColor',fcolor,'YMinorTick','on','XMinorTick','on','fontsize',26,'LineWidth',2)
    set(xl, 'fontweight','bold','color',fcolor,'fontsize',32);
    set(yl,'fontweight','bold','color',fcolor,'fontsize',32); set(zl,'fontweight','bold','color',fcolor,'fontsize',32);
    hold on
    sizeY = size(image.RawStack,2)*raw.dataKey(9,1);
    
    
    boxC = [.5 .7 .9];
    
    fill3(RXtb,RYtb,RZsb,boxC,'FaceAlpha',.1,'EdgeAlpha',1,'LineStyle','--','EdgeColor',boxC,'LineWidth',2)
    fill3(RXtb,RYtb,RZst,boxC,'FaceAlpha',.1,'EdgeAlpha',1,'LineStyle','--','EdgeColor',boxC,'LineWidth',2)
    fill3(RXfb,RYf,RZss,boxC,'FaceAlpha',.1,'EdgeAlpha',1,'LineStyle','--','EdgeColor',boxC,'LineWidth',2)
    fill3(RXfb,RYb,RZss,boxC,'FaceAlpha',.1,'EdgeAlpha',1,'LineStyle','--','EdgeColor',boxC,'LineWidth',2)
    fill3(RXs1,RYss,RZss,boxC,'FaceAlpha',.1,'EdgeAlpha',1,'LineStyle','--','EdgeColor',boxC,'LineWidth',2)
    fill3(RXs2,RYss,RZss,boxC,'FaceAlpha',.1,'EdgeAlpha',1,'LineStyle','--','EdgeColor',boxC,'LineWidth',2)
    
    fill3([60 60 100 100], [35 35 30 30], [bot top top bot],[0 1 0],'FaceAlpha',.3,'EdgeAlpha',1,'LineStyle','-','EdgeColor','green','LineWidth',3)
    
    xlim([figLimits(1,1)-2 figLimits(1,2)+2])
    ylim([figLimits(1,3) figLimits(1,4) ])
    zlim([3 12])
    
    savefile = [strcat(filePath,'\MethodsExpFigs\') fcolor ' on ' bcolor 'fig9 small.tif'];
    export_fig(f9,savefile,'-native');
    
    %%
    %Figure 10b - Quiver Plot Full With Profile plane Trace Viewed from XZ
    
    az= -0;
    ele= 0;
    
    
    
    f9 = figure;
    set(gcf,'position',[100,100,1200,600])
    hold on
            cm(1:3,1) = ([0.8500    0.3250    0.0980]);
        cm(1:3,2) = ([0.9290    0.6940    0.1250]);
    for i = [1 2]%:size(plane.final,2)

        quiver3(m3.refSC(plane.final(1:nnz(plane.final(:,i)),i),1),(size(image.Area,1)*raw.dataKey(9,1))-m3.refSC(plane.final(1:nnz(plane.final(:,i)),i),2),m3.refSC(plane.final(1:nnz(plane.final(:,i)),i),3),m3.disp(plane.final(1:nnz(plane.final(:,i)),i),1),m3.disp(plane.final(1:nnz(plane.final(:,i)),i),2)*-1,m3.disp(plane.final(1:nnz(plane.final(:,i)),i),3),0,'Color',cm(1:3,i))
        bots(i,1) = min(r.r(r.r(plane.final(1:nnz(plane.final(:,i)),i),3)>0,3));
    end
    xl = xlabel('X (\mum)');
    yl = ylabel('Y (\mum)');
    zl = zlabel('Z (\mum)');
    view([az ele])
    set(gca,'Color',bcolor)
    set(gca,'XColor',fcolor,'YColor',fcolor,'ZColor',fcolor,'YMinorTick','on','XMinorTick','on','fontsize',26,'LineWidth',2)
    set(xl, 'fontweight','bold','color',fcolor,'fontsize',32);
    set(yl,'fontweight','bold','color',fcolor,'fontsize',32); set(zl,'fontweight','bold','color',fcolor,'fontsize',32);
    hold on
    sizeY = size(image.RawStack,2)*raw.dataKey(9,1);
    RXtb = [figLimits(1,1) figLimits(1,2) figLimits(1,2) figLimits(1,1) figLimits(1,1)];
    RXfb = [figLimits(1,1) figLimits(1,1) figLimits(1,2) figLimits(1,2) figLimits(1,1)];
    RXs1 = [figLimits(1,1) figLimits(1,1) figLimits(1,1) figLimits(1,1) figLimits(1,1)];
    RXs2 = [figLimits(1,2) figLimits(1,2) figLimits(1,2) figLimits(1,2) figLimits(1,2)];
    
    RYtb = [figLimits(1,3) figLimits(1,3) figLimits(1,4) figLimits(1,4) figLimits(1,3)];
    RYf = [figLimits(1,3) figLimits(1,3) figLimits(1,3) figLimits(1,3) figLimits(1,3)];
    RYb = [figLimits(1,4) figLimits(1,4) figLimits(1,4) figLimits(1,4) figLimits(1,4)];
    RYss = [figLimits(1,3) figLimits(1,3) figLimits(1,4) figLimits(1,4) figLimits(1,3)];
    
    
    
    bot = min(bots)
    top = max(r.r(:,3))
    RZsb = [bot bot bot bot bot];
    RZst = [top top top top top];
    RZss = [bot top top bot bot];
    
    boxC = [.5 .7 .9];
    
    %plot3(RXtb,RYtb,RZsb,'LineStyle','--','Color',boxC,'LineWidth',2)
    %plot3(RXtb,RYtb,RZst,'LineStyle','--','Color',boxC,'LineWidth',2)
    plot3(RXfb,RYf,RZss,'LineStyle','--','Color',boxC,'LineWidth',2)
    plot3(RXfb,RYb,RZss,'LineStyle','--','Color',boxC,'LineWidth',2)
    %plot3(RXs1,RYss,RZss,'LineStyle','--','Color',boxC,'LineWidth',2)
    %plot3(RXs2,RYss,RZss,'LineStyle','--','Color',boxC,'LineWidth',2)
    
    plot3([60 60 100 100 60], [35 35 30 30 35], [bot top top bot bot],'Color','green','LineStyle','-','LineWidth',3)
    
    savefile = [strcat(filePath,'\MethodsExpFigs\') fcolor ' on ' bcolor 'fig10 Side.tif'];
    export_fig(f9,savefile,'-native');
    %%
    %Figure 10c - Quiver Plot Full With Profile plane Trace Viewed from XY
    
    az= 0;
    ele= 90;
    
    
    
    f9 = figure;
    set(gcf,'position',[0,0,500,500])
    hold on
    plot3(RXtb,RYtb,RZsb,'LineStyle','--','Color',boxC,'LineWidth',2)
    plot3(RXtb,RYtb,RZst,'LineStyle','--','Color',boxC,'LineWidth',2)
    plot3(RXfb,RYf,RZss,'LineStyle','--','Color',boxC,'LineWidth',2)
    plot3(RXfb,RYb,RZss,'LineStyle','--','Color',boxC,'LineWidth',2)
    plot3(RXs1,RYss,RZss,'LineStyle','--','Color',boxC,'LineWidth',2)
    plot3(RXs2,RYss,RZss,'LineStyle','--','Color',boxC,'LineWidth',2)
    
    plot3([60 60 100 100 60], [35 35 30 30 35], [bot top top bot bot],'Color','green','LineStyle','-','LineWidth',3)
    for i = [2]%:size(plane.final,2)
        quiver3(m3.refSC(plane.final(1:nnz(plane.final(:,i)),i),1),(size(image.Area,1)*raw.dataKey(9,1))-m3.refSC(plane.final(1:nnz(plane.final(:,i)),i),2),m3.refSC(plane.final(1:nnz(plane.final(:,i)),i),3),m3.disp(plane.final(1:nnz(plane.final(:,i)),i),1),m3.disp(plane.final(1:nnz(plane.final(:,i)),i),2)*-1,m3.disp(plane.final(1:nnz(plane.final(:,i)),i),3),2,'Color',[0.9290    0.6940    0.1250])
        bots(i,1) = min(r.r(r.r(plane.final(1:nnz(plane.final(:,i)),i),3)>0,3));
    end
    xl = xlabel('X (\mum)');
    yl = ylabel('Y (\mum)');
    zl = zlabel('Z (\mum)');
    view([az ele])
    set(gca,'Color',bcolor)
    set(gca,'XColor',fcolor,'YColor',fcolor,'ZColor',fcolor,'YMinorTick','on','XMinorTick','on','fontsize',26,'LineWidth',2)
    set(xl, 'fontweight','bold','color',fcolor,'fontsize',32);
    set(yl,'fontweight','bold','color',fcolor,'fontsize',32); set(zl,'fontweight','bold','color',fcolor,'fontsize',32);
    hold on
    sizeY = size(image.RawStack,2)*raw.dataKey(9,1);

    
    boxC = [.5 .7 .9];
    
    ylim([0 80])
    

    
    savefile = [strcat(filePath,'\MethodsExpFigs\') fcolor ' on ' bcolor 'fig10 top.tif'];
    export_fig(f9,savefile,'-native');
    
    
    
end

