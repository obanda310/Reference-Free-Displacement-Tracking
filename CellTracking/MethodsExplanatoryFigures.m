%% Methods script for creating figures to accompany discriptions of reference generation and displacement calculations
load('3DnormalData.mat')
mkdir('MethodsExpFigs')
%%
close all
[~,figBounds] =  imcrop(HeatMap3(:,:,:,1));
figLimits(1,1) = figBounds(1,1) *xyScale;
figLimits(1,2) = (figBounds(1,1) + figBounds(1,3))*xyScale;
figLimits(1,3) = figBounds(1,2) *xyScale;
figLimits(1,4) = (figBounds(1,2) + figBounds(1,4))*xyScale;
% Figure 1 - Forming columns and references for XY data: initial references
%%
save('figBounds.mat','figLimits','figBounds')
res2 = permute(res,[2 1 3]);
res3 = res2(figBounds(1,1):figBounds(1,1)+figBounds(1,3),figBounds(1,2):figBounds(1,2)+figBounds(1,4),:);
%ShowStack(res3)
StackFile = [filePath,'\MethodsExpFigs\','CroppedRegion.tif'];
for i = 1:size(res3,3)
imwrite(uint16(res3(:,:,i)),StackFile,'WriteMode','append');
end
%%
az = -110;
ele = 50;
%%
f0 = figure;

hold on
for i = 1:size(book1,3)    
    %plot3(ones(size(book1,2),1)*book1(1,1,i),size(res,2)*xyScale-ones(size(book1,2),1)*book1(2,1,i),((1:1:size(book1,2))*zScale)','r')
    scatter3(book1(1,:,i),size(res,2)*xyScale-book1(2,:,i),(1:1:size(book1,2))*zScale,'.')
end
xlim([figLimits(1,1) figLimits(1,2)])
ylim([figLimits(1,3) figLimits(1,4)])
xlabel('X (\mum)')
ylabel('Y (\mum)')
zlabel('Z (\mum)')

view([az ele])

    savefile = [strcat(filePath,'\MethodsExpFigs\') 'fig0.tif'];
    export_fig(f0,savefile,'-native');
    %%
f0b = figure;

hold on
for i = 1:size(book1,3)    
    %plot3(ones(size(book1,2),1)*book1(1,1,i),size(res,2)*xyScale-ones(size(book1,2),1)*book1(2,1,i),((1:1:size(book1,2))*zScale)','r')
    scatter3(book1(1,:,i),size(res,2)*xyScale-book1(2,:,i),(1:1:size(book1,2))*zScale,'.','b')
end
xlim([figLimits(1,1) figLimits(1,2)])
ylim([figLimits(1,3) figLimits(1,4)])
xlabel('X (\mum)')
ylabel('Y (\mum)')
zlabel('Z (\mum)')

view([az ele])

    savefile = [strcat(filePath,'\MethodsExpFigs\') 'fig0b.tif'];
    export_fig(f0b,savefile,'-native');
    
%%

f1 = figure;

hold on
for i = 1:size(book1,3)    
    plot3(ones(size(book1,2),1)*book1(1,1,i),size(res,2)*xyScale-ones(size(book1,2),1)*book1(2,1,i),((1:1:size(book1,2))*zScale)','r')
    scatter3(book1(1,:,i),size(res,2)*xyScale-book1(2,:,i),(1:1:size(book1,2))*zScale,'.')
end
xlim([figLimits(1,1) figLimits(1,2)])
ylim([figLimits(1,3) figLimits(1,4)])
xlabel('X (\mum)')
ylabel('Y (\mum)')
zlabel('Z (\mum)')

view([az ele])

    savefile = [strcat(filePath,'\MethodsExpFigs\') 'fig1.tif'];
    export_fig(f1,savefile,'-native');
%%
book1(22,:,:) = book1(1,:,:) - book1(13,:,:)*xyScale;
book1(23,:,:) = book1(2,:,:) - book1(14,:,:)*xyScale;
%Figure 2 - Tilt Corrected
f2 = figure;
hold on
for i = 1:size(book1,3)    
    plot3(book1(22,:,i),size(res,2)*xyScale-book1(23,:,i),((1:1:size(book1,2))*zScale)','r')
    scatter3(book1(1,:,i),size(res,2)*xyScale-book1(2,:,i),(1:1:size(book1,2))*zScale,'.')
end
xlim([figLimits(1,1) figLimits(1,2)])
ylim([figLimits(1,3) figLimits(1,4)])
xlabel('X (\mum)')
ylabel('Y (\mum)')
zlabel('Z (\mum)')
view([az ele])
savefile = [strcat(filePath,'\MethodsExpFigs\') 'fig2.tif'];
    export_fig(f2,savefile,'-native');
%%
%Figure 3 - plus Quiver

f3 = figure;
hold on
for i = 1:size(book1,3)    
    plot3(book1(22,:,i),size(res,2)*xyScale-book1(23,:,i),((1:1:size(book1,2))*zScale)','r')
    scatter3(book1(1,:,i),size(res,2)*xyScale-book1(2,:,i),(1:1:size(book1,2))*zScale,'.')
    quiver3(book1(22,:,i),size(res,2)*xyScale-book1(23,:,i),((1:1:size(book1,2))*zScale),book1(13,:,i)*xyScale,-1*book1(14,:,i)*xyScale,zeros(size(book1,2),1)',0)
end

xlim([figLimits(1,1) figLimits(1,2)])
ylim([figLimits(1,3) figLimits(1,4)])
xlabel('X (\mum)')
ylabel('Y (\mum)')
zlabel('Z (\mum)')
view([az ele])
savefile = [strcat(filePath,'\MethodsExpFigs\') 'fig3.tif'];
    export_fig(f3,savefile,'-native');
%%
%Figure 4 - Quiver plot only

f4 = figure;
hold on
for i = 1:size(book1,3)    
    quiver3(book1(22,:,i),size(res,2)*xyScale-book1(23,:,i),((1:1:size(book1,2))*zScale),book1(13,:,i)*xyScale,-1*book1(14,:,i)*xyScale,zeros(size(book1,2),1)',0)
end

xlim([figLimits(1,1) figLimits(1,2)])
ylim([figLimits(1,3) figLimits(1,4)])
xlabel('X (\mum)')
ylabel('Y (\mum)')
zlabel('Z (\mum)')
view([az ele])
savefile = [strcat(filePath,'\MethodsExpFigs\') 'fig4.tif'];
    export_fig(f4,savefile,'-native');

%%
%Change elevation to better view 3D detections
az = -110;
ele = 15;

%% Gather detection rows in window
rZone = unique(r(r(:,1)>figLimits(1,1)&r(:,1)<figLimits(1,2)&r(:,2)<size(res,2)*xyScale-figLimits(1,3)&r(:,2)>size(res,2)*xyScale-figLimits(1,4),8));
rZonePts = find(r(:,1)>figLimits(1,1)&r(:,1)<figLimits(1,2)&r(:,2)<size(res,2)*xyScale-figLimits(1,3)&r(:,2)>size(res,2)*xyScale-figLimits(1,4));
%%
%Figure 5 - 3D Detections

f5 = figure;
hold on
for i = 1:size(book1,3)    
    
    scatter3(book1(1,:,i),size(res,2)*xyScale-book1(2,:,i),(1:1:size(book1,2))*zScale,'.')
end

for i =1:size(planesFinal,2)
scatter3(r(planesFinal(1:nnz(planesFinal(:,i)),i),1),size(res,2)*xyScale-r(planesFinal(1:nnz(planesFinal(:,i)),i),2),r(planesFinal(1:nnz(planesFinal(:,i)),i),3))
end
xlim([figLimits(1,1) figLimits(1,2)])
ylim([figLimits(1,3) figLimits(1,4)])
xlabel('X (\mum)')
ylabel('Y (\mum)')
zlabel('Z (\mum)')
view([az ele])


    savefile = [strcat(filePath,'\MethodsExpFigs\') 'fig5.tif'];
    export_fig(f5,savefile,'-native');
%%
%Figure 6a - 3D Detections and Row Assignment

f6a = figure;
hold on
colorMapFig7 = brewermap(size(rZone,1),'*spectral');
clear temp
for i = 1:size(rZone,1)
    scatter3(r(rows(rZone(i,1),1:nnz(rows(rZone(i,1),:))),1),(size(imageBinary,1)*xyScale)-r(rows(rZone(i,1),1:nnz(rows(rZone(i,1),:))),2),r(rows(rZone(i,1),1:nnz(rows(rZone(i,1),:))),3),'MarkerEdgeColor',[colorMapFig7(i,1:3)])
    temp(1:nnz(rows(rZone(i,1),:)),1:3,i) = sortrows(r(rows(rZone(i,1),1:nnz(rows(rZone(i,1),:))),1:3),1);
    plot3(temp(:,1,i),(size(imageBinary,1)*xyScale)-temp(:,2,i),temp(:,3,i),'-','Color',[colorMapFig7(i,1:3)])
    %plot3([rowFits3(i,1,1) rowFits3(i,1,2)],[(size(imageBinary,1)*xyScale)-rowFits3(i,2,1) (size(imageBinary,1)*xyScale)-rowFits3(i,2,2)],[rowFits3(i,3,1) rowFits3(i,3,2)],'Color',[colorMapFig7(i,1:3)])
end
scatter3(0,0,0)

xlim([figLimits(1,1) figLimits(1,2)])
ylim([figLimits(1,3) figLimits(1,4)])
xlabel('X (\mum)')
ylabel('Y (\mum)')
zlabel('Z (\mum)')
view([az ele])
temp(temp==0) = NaN;

    savefile = [strcat(filePath,'\MethodsExpFigs\') 'fig6a.tif'];
    export_fig(f6a,savefile,'-native');
    %%
%Figure 6b - 3D Detections and Row Fits

az2=-130;
ele2=30;

f6b = figure;
hold on
colorMapFig7 = brewermap(size(rZone,1),'*spectral');
rowsFilt = rows;
for i = 1:size(rZone,1)
    rowsFilt(rZone(i,1),:) = 0;
    rowsFilt(rZone(i,1),1:nnz(intersect(rZonePts,rows(rZone(i,1),:)))) = intersect(rZonePts,rows(rZone(i,1),:));
end

for i = 1:size(rZone,1)
    scatter3(r(rows(rZone(i,1),1:nnz(rows(rZone(i,1),:))),1),(size(imageBinary,1)*xyScale)-r(rows(rZone(i,1),1:nnz(rows(rZone(i,1),:))),2),r(rows(rZone(i,1),1:nnz(rows(rZone(i,1),:))),3),'MarkerEdgeColor',[.5 .5 .5])
    plot3(temp(:,1,i),(size(imageBinary,1)*xyScale)-temp(:,2,i),temp(:,3,i),'-','Color',[.5 .5 .5])
    plot3([rowFits3(rZone(i,1),1,1) rowFits3(rZone(i,1),1,2)],[(size(imageBinary,1)*xyScale)-rowFits3(rZone(i,1),2,1) (size(imageBinary,1)*xyScale)-rowFits3(rZone(i,1),2,2)],[rowFits3(rZone(i,1),3,1) rowFits3(rZone(i,1),3,2)],'r')
end

for i = 1:size(rZone,1)
    scatter3(r(rowsFilt(rZone(i,1),1:nnz(rowsFilt(rZone(i,1),:))),1),(size(imageBinary,1)*xyScale)-r(rowsFilt(rZone(i,1),1:nnz(rowsFilt(rZone(i,1),:))),2),r(rowsFilt(rZone(i,1),1:nnz(rowsFilt(rZone(i,1),:))),3),'MarkerEdgeColor',[colorMapFig7(i,1:3)])
    temp2(:,1:3) = sortrows(r(rowsFilt(rZone(i,1),1:nnz(rowsFilt(rZone(i,1),:))),1:3),1);
    plot3(temp2(:,1),(size(imageBinary,1)*xyScale)-temp2(:,2),temp2(:,3),'Color',[colorMapFig7(i,1:3)])
end

%xlim([figLimits(1,1) figLimits(1,2)])
%ylim([figLimits(1,3) figLimits(1,4)])
xlabel('X (\mum)')
ylabel('Y (\mum)')
zlabel('Z (\mum)')
view([az2 ele2])


    savefile = [strcat(filePath,'\MethodsExpFigs\') 'fig6b.tif'];
    export_fig(f6b,savefile,'-native');
        %%
%Figure 6C - 3D Detections and Row Fits

f6c = figure;
hold on
colorMapFig7 = brewermap(size(rZone,1),'*spectral');
rowsFilt = rows;
for i = 1:size(rZone,1)
    rowsFilt(rZone(i,1),:) = 0;
    rowsFilt(rZone(i,1),1:nnz(intersect(rZonePts,rows(rZone(i,1),:)))) = intersect(rZonePts,rows(rZone(i,1),:));
end

for i = 1:size(rZone,1)
    scatter3(r(rows(rZone(i,1),1:nnz(rows(rZone(i,1),:))),1),(size(imageBinary,1)*xyScale)-r(rows(rZone(i,1),1:nnz(rows(rZone(i,1),:))),2),r(rows(rZone(i,1),1:nnz(rows(rZone(i,1),:))),3),'MarkerEdgeColor',[.5 .5 .5])
    plot3(temp(:,1,i),(size(imageBinary,1)*xyScale)-temp(:,2,i),temp(:,3,i),'-','Color',[.5 .5 .5])
    plot3([rowFits3(rZone(i,1),1,1) rowFits3(rZone(i,1),1,2)],[(size(imageBinary,1)*xyScale)-rowFits3(rZone(i,1),2,1) (size(imageBinary,1)*xyScale)-rowFits3(rZone(i,1),2,2)],[rowFits3(rZone(i,1),3,1) rowFits3(rZone(i,1),3,2)],'r')
end

for i = 1:size(rZone,1)
    scatter3(r(rowsFilt(rZone(i,1),1:nnz(rowsFilt(rZone(i,1),:))),1),(size(imageBinary,1)*xyScale)-r(rowsFilt(rZone(i,1),1:nnz(rowsFilt(rZone(i,1),:))),2),r(rowsFilt(rZone(i,1),1:nnz(rowsFilt(rZone(i,1),:))),3),'MarkerEdgeColor',[colorMapFig7(i,1:3)])
    temp2(:,1:3) = sortrows(r(rowsFilt(rZone(i,1),1:nnz(rowsFilt(rZone(i,1),:))),1:3),1);
    plot3(temp2(:,1),(size(imageBinary,1)*xyScale)-temp2(:,2),temp2(:,3),'Color',[colorMapFig7(i,1:3)])
end

%xlim([figLimits(1,1) figLimits(1,2)])
%ylim([figLimits(1,3) figLimits(1,4)])
xlim([figLimits(1,1) figLimits(1,2)])
ylim([figLimits(1,3) figLimits(1,4)])
xlabel('X (\mum)')
ylabel('Y (\mum)')
zlabel('Z (\mum)')
view([az ele])


    savefile = [strcat(filePath,'\MethodsExpFigs\') 'fig6c.tif'];
    export_fig(f6c,savefile,'-native');
%%
%Figure 7a - Row Fits/Column Fits
f7a = figure;
hold on
colorMapFig7 = brewermap(size(rows,1),'*spectral');

for i = 1:size(book1,3)    
    plot3(book1(22,:,i),size(res,2)*xyScale-book1(23,:,i),((1:1:size(book1,2))*zScale)','b')
end

for i = 1:size(rowFits3,1)
    
    plot3([rowFits3(i,1,1) rowFits3(i,1,2)],[(size(imageBinary,1)*xyScale)-rowFits3(i,2,1) (size(imageBinary,1)*xyScale)-rowFits3(i,2,2)],[rowFits3(i,3,1) rowFits3(i,3,2)],'Color','r')
end
scatter3(rRef3b(:,1),(size(imageBinary,1)*xyScale)-rRef3b(:,2),rRef3b(:,3),'x','r')
scatter3(0,0,0)

xlim([figLimits(1,1) figLimits(1,2)])
ylim([figLimits(1,3) figLimits(1,4)])
xlabel('X (\mum)')
ylabel('Y (\mum)')
zlabel('Z (\mum)')
view([az ele])

    savefile = [strcat(filePath,'\MethodsExpFigs\') 'fig7a.tif'];
    export_fig(f7a,savefile,'-native');
%%
%Figure 7b - 3D Detections/Row Fits/Column Fits

f7b = figure;
hold on
colorMapFig7 = brewermap(size(rows,1),'*spectral');

for i = 1:size(rowFits3,1)
    scatter3(r(rows(i,1:nnz(rows(i,:))),1),(size(imageBinary,1)*xyScale)-r(rows(i,1:nnz(rows(i,:))),2),r(rows(i,1:nnz(rows(i,:))),3),'MarkerEdgeColor',[colorMapFig7(i,1:3)])
    plot3([rowFits3(i,1,1) rowFits3(i,1,2)],[(size(imageBinary,1)*xyScale)-rowFits3(i,2,1) (size(imageBinary,1)*xyScale)-rowFits3(i,2,2)],[rowFits3(i,3,1) rowFits3(i,3,2)],'Color',[colorMapFig7(i,1:3)])
end
scatter3(rRef3b(:,1),(size(imageBinary,1)*xyScale)-rRef3b(:,2),rRef3b(:,3),'x','r')
scatter3(0,0,0)

xlim([figLimits(1,1) figLimits(1,2)])
ylim([figLimits(1,3) figLimits(1,4)])
xlabel('X (\mum)')
ylabel('Y (\mum)')
zlabel('Z (\mum)')
view([az ele])

    savefile = [strcat(filePath,'\MethodsExpFigs\') 'fig7b.tif'];
    export_fig(f7b,savefile,'-native');
%%
%Figure 8 - 3D Detections/Row Fits/Quiver
f8 = figure;
hold on
colorMapFig7 = brewermap(size(rows,1),'*spectral');

for i = 1:size(rowFits3,1)
    scatter3(r(rows(i,1:nnz(rows(i,:))),1),(size(imageBinary,1)*xyScale)-r(rows(i,1:nnz(rows(i,:))),2),r(rows(i,1:nnz(rows(i,:))),3),'MarkerEdgeColor',[colorMapFig7(i,1:3)])
    plot3([rowFits3(i,1,1) rowFits3(i,1,2)],[(size(imageBinary,1)*xyScale)-rowFits3(i,2,1) (size(imageBinary,1)*xyScale)-rowFits3(i,2,2)],[rowFits3(i,3,1) rowFits3(i,3,2)],'Color',[colorMapFig7(i,1:3)])
end
scatter3(0,0,0)
quiver3(rRef3b(:,1),(size(imageBinary,1)*xyScale)-rRef3b(:,2),rRef3b(:,3),rDisp3(:,1),rDisp3(:,2)*-1,rDisp3(:,3),0)
xlim([figLimits(1,1) figLimits(1,2)])
ylim([figLimits(1,3) figLimits(1,4)])
xlabel('X (\mum)')
ylabel('Y (\mum)')
zlabel('Z (\mum)')
view([az ele])


    savefile = [strcat(filePath,'\MethodsExpFigs\') 'fig8.tif'];
    export_fig(f8,savefile,'-native');
    
%%
%Figure 9 - Quiver Plot Full
az= -30;
ele= 30;

f9 = figure;
hold on
for i = 1:size(planesFinal,2)
quiver3(rRef3b(planesFinal(1:nnz(planesFinal(:,i)),i),1),(size(imageBinary,1)*xyScale)-rRef3b(planesFinal(1:nnz(planesFinal(:,i)),i),2),rRef3b(planesFinal(1:nnz(planesFinal(:,i)),i),3),rDisp3(planesFinal(1:nnz(planesFinal(:,i)),i),1),rDisp3(planesFinal(1:nnz(planesFinal(:,i)),i),2)*-1,rDisp3(planesFinal(1:nnz(planesFinal(:,i)),i),3),0)
end
xlabel('X (\mum)')
ylabel('Y (\mum)')
zlabel('Z (\mum)')
    view([az ele])
    
     savefile = [strcat(filePath,'\MethodsExpFigs\') 'fig9.tif'];
    export_fig(f9,savefile,'-native');