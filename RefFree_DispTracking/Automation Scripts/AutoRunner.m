%% Automated Script
clear all
close all
[ListName,ListPath] = uigetfile;
fullListName = strcat(ListPath,ListName);
load(strcat(ListPath,ListName),'dirList')
parts = strsplit(ListName, 'List');
prefix = parts{1};
set(0,'defaultfigurecolor',[1 1 1])
%% For updating paths (when files are moved)
%set a as new path root
%set b as old path root
ak = strfind(fullListName,ListName);

a = [fullListName(1:ak-1) 'Group'];

Schar = char(dirList{:});
all_rows_same = all(diff(Schar, [], 1) == 0, 1);
common_cols = find(~all_rows_same, 1, 'first');
if isempty(common_cols)
    b = '?'
else
    b = dirList{1}(1:common_cols-1);
end

for i = 1:size(dirList,1)
    
    dirList{i,1} = strrep(dirList{i,1},b,a);
end
%%
save(fullListName,'dirList')

%% For clearing old files
for j = 1:2
    for i = 1:size(dirList,1)
        cleanOld(dirList{i,1})
        disp(num2str(i))
    end
end
%% For Running Tracking Script
for i = 1:size(dirList,1)
    try
        tracking(dirList{i,1})
        disp(num2str(i))
    catch
        disp('Tracking Failed')
        i
    end
end
%% For Running XY Script
for i = 1:size(dirList,1)
    try
        dispShear(dirList{i,1})
        disp(num2str(i))
    catch
        disp('Shear Failed')
        i
    end
end

%% For Running 3D Script
for i = 1:size(dirList,1)
    try
        disp3D(dirList{i,1})
        disp(num2str(i))
    catch
        disp('3D Failed')
        i
    end
end
%% For Running Spherical Indentation Script
for i = 1:size(dirList,1)
    try
        SphericalIndent(dirList{i,1})
        disp(num2str(i))
    catch
        disp('Indentation Failed')
        i
    end
end

%% For calculating mean deviations in a dataset
clear StdVals
for i2 = 1:size(dirList,1)
    cd(dirList{i2,1})
    load('3Ddata.mat')
    xVals = shear.ltdX(:,shear.noCellTraj);
    xVals(xVals == 0) = NaN;
    
    yVals = shear.ltdY(:,shear.noCellTraj);
    yVals(yVals == 0) = NaN;
    
    xyVals = squeeze(shear.ltdXY(:,shear.noCellTraj));
    xyVals = cat(1,xyVals,-1*xyVals);
    xyVals(xyVals == 0) = NaN;
    
    StdVals(i2,1) = std(xVals(:),'omitnan')*1000;
    StdVals(i2,2) = std(yVals(:),'omitnan')*1000;
    StdVals(i2,3) = std(xyVals(:),'omitnan')*1000;
    StdVals(i2,4) = m3.noiseCutoff*1000;
    
    i2
    
end
%%
xyzStats(1,1) = mean(StdVals(:,3)*2);
xyzStats(1,2) = std(StdVals(:,3)*2);

xyzStats(1,3) = mean(StdVals((StdVals(:,3)<(xyzStats(1,1)+2*xyzStats(1,2))),3)*2);
xyzStats(1,4) = std(StdVals((StdVals(:,3)<(xyzStats(1,1)+2*xyzStats(1,2))),3)*2);

xyzStats(2,1) = mean(StdVals(:,4));
xyzStats(2,2) = std(StdVals(:,4));

xyzStats(2,3) = mean(StdVals((StdVals(:,4)<(xyzStats(2,1)+2*xyzStats(2,2))),4));
xyzStats(2,4) = std(StdVals((StdVals(:,4)<(xyzStats(2,1)+2*xyzStats(2,2))),4));


%% For Running Z Script
for i = 1:size(dirList,1)
    cd(dirList{i,1})
    %First check to see if the stack to analyze has been stored previously
    clear files
    files = dir('*.mat');
    if size(files,1)>=1
        for k = 1:length(files)
            current=files(k).name;
            if size(current,2)>12
                check(k)=strcmp(current(end-12:end),'StackName.mat');
            else
                check(k) = 0;
            end
        end
    end
    
    %If the Stack has not been assigned previously, manually assign it
    if size(find(check),2)==0
        [StackName,StackPath] = uigetfile('*.tif');
        save('StackName.mat','StackName')
    else
        load('StackName.mat','StackName')
        if StackName == 0
            [StackName,StackPath] = uigetfile('*.tif');
            save('StackName.mat','StackName')
        end
    end
    
end
%% Make sure every folder has a rowV file
[rowVName,rowVPath] = uigetfile;
load(strcat(rowVPath,rowVName),'rowV')
for i = 1:size(dirList,1)
    cd(dirList{i,1})
    save('rowV.mat','rowV')
end

%% For Running Z Script
for i = 1:size(dirList,1)
    try
        disp(dirList{i,1})
        disp3D(dirList{i,1})
        disp(dirList{i,1})
        i
    catch
        disp('Normal Failed')
        i
    end
end

%% For Running Z Script Outputs Only (faster alternative)
for i = 1:size(dirList,1)
    disp3DOutputs(dirList{i,1})
    i
end


%%
for i = 1:size(dirList,1)
    disp(dirList{i,1})
    cd(dirList{i,1});
    load('3Ddata.mat')
    filePath = cd;
    folderName = 'Profile Data';
    mkdir(filePath,folderName)
    save('Profile Data\vqZ.mat','vqN','vq3','image','HeatMapN','HeatMap3')
end

%% Collect vqXY/vqZ totals for shear v normal graphs
%goodData = [1;2;3;5;6;14;17;23;31;34;43;44;47;49;50;48;11;12;19;25;28;36];
%goodData = [1:1:51]';

goodData = [1:1:14]';
clear vqT
for i2 = 1:size(goodData,1)
    try
        clear planesLoc3
        disp(dirList{goodData(i2,1),1})
        cd(dirList{goodData(i2,1),1});
        load('3Ddata.mat')
        load('Profile Data\vqXY.mat');
        load('Profile Data\HeatMapXY.mat');
        zTarget = 7;
        for j = 1:size(planesGroups,1)
            planesLoc3(j) = mean(planesLoc2(1,planesGroups(j,1:nnz(planesGroups(j,:)))));
        end
        zPlane = find(abs(planesLoc3-zTarget) == min(abs(planesLoc3-zTarget)),1,'first');
        vqXY(isnan(vqXY)) = 0;
        vqXY = imresize(vqXY,[size(vq3,1) size(vq3,2)]);
        i2
        vqT(i2,1) = sum(sum(vqXY));
        vq3(isnan(vq3)) = 0;
        vqT(i2,2) = sum(sum(abs(vq3(:,:,zPlane))));
        vqT(i2,3) = zPlane;
        vqT(i2,4) = planesLoc2(zPlane);
        vqT(i2,5) = sum(sum(image.Area ==0))*raw.dataKey(9,1)^2;
        vqT(i2,6) = max(max(vqXY));
        vqT(i2,7) = max(max(abs(vq3(:,:,zPlane))));
    catch
        
    end
end
%%
vqT(:,8) = vqT(:,1)./vqT(:,2);
vqT(:,9) = 1; vqT(:,1)>prctile(vqT(:,1),0);
vqTM = mean(vqT(vqT(:,9)==1,8));
vqTS = std(vqT(vqT(:,9)==1,8));
save([ListPath prefix 'ShearNormalStats.mat'],'vqT','vqTS','vqTM')

%%
figure
scatter(vqT(:,1),vqT(:,2))
figure
scatter(vqT(:,5),vqT(:,1))
figure
scatter(vqT(:,5),vqT(:,2))
figure
scatter(vqT(:,6),vqT(:,7))
%% Imbalance Calculator
clear ftotals
for i = 1:size(dirList,1)
    clear totals
    totals = imbalance(dirList{i});
    ftotals(i,:) = totals;
    i
end
%%
ftotals(isnan(ftotals))=0;
(sum(ftotals(:,2)) + sum(ftotals(:,4)))/sum(ftotals(:,6))
%%

Ys(2,:) = ftotals(:,2)+ftotals(:,4);
Ys(3,:) = ftotals(:,1)+ftotals(:,3);
BarGrapher(Ys)

%% For Running Profile Script
clear caught goodData
goodData = [1;2;3;5;6;14;17;23;31;34;43;44;47;49;50;48;11;12;19;25;28;36];
for i = 1:size(goodData,1)
    try
        [normXY,normZ,normAxis,cell_boundary] = profFunc(dirList{goodData(i,1),1});
        profBook(1,1:size(normAxis,2),i) = normAxis;
        profBook(2,1:size(normXY,2),i) = normXY;
        profBook(3,1:size(normZ,2),i) = normZ;
        cb(i) = cell_boundary;
        save('Profile Data\ProfileOutputData.mat','normXY','normZ','normAxis')
    catch
        caught(i) = 1;
        i
        dirList{goodData(i,:),1}
    end
end
cbf = mean(cb);
save(strcat(ListPath,prefix,'Profiles.mat'),'profBook','cb','cbf')

%%
prefix = '\Gels';
cb(cb<25 | cb>75)=NaN;
cbf2 = mean(cb,'omitnan')/100;
aScale = 1/cbf2;
clear profBookCat
profBookCat = profBook;
clear keep keep2 keep3 map profBook2 profBook3
for i = 1:size(profBookCat,3)
    %max(profBookCat(2,:,i))
    if max(profBookCat(2,:,i))>1
        keep(1,i) = 1;
        keep(2,i) = max(profBookCat(2,:,i));
        profBook2(1:size(profBookCat(:,:,i),1),1:size(profBookCat(:,:,i),2),i) = profBookCat(:,:,i);
    end
end
profBook2(:,:,~keep(1,:)) = [];
keep2 = keep(2,:);
keep2(:,~keep(1,:)) = [];
[keep3,sortIdx] = sort(keep2);

colOptions{1,1} = 'white';
colOptions{2,1} = 'black';
colOptions{1,2} = 'black';
colOptions{2,2} = 'white';
map = brewermap(size(profBook2,3),'*spectral');

for i = 1:size(colOptions,2)
    fcolor = colOptions{1,i};
    bcolor = colOptions{2,i};
    
    set(0,'defaultfigurecolor',bcolor)
    ProfileOverlays = figure;
    hold on
    for i = 1:size(profBook2,3)
        plot(profBook2(1,:,sortIdx(i))*aScale,profBook2(2,:,sortIdx(i)),'color',[map(i,1:3)])
        plot(profBook2(1,:,sortIdx(i))*aScale,profBook2(3,:,sortIdx(i)),'color',[map(i,1:3)],'linestyle','-.')
    end
    p1= plot(profBook2(1,:,sortIdx(i))*aScale,mean(profBook2(2,:,sortIdx(:)),3,'omitnan'),'color',fcolor,'linewidth',3);
    p2=plot(profBook2(1,:,sortIdx(i))*aScale,mean(profBook2(3,:,sortIdx(:)),3,'omitnan'),'color',fcolor,'linestyle','-.','linewidth',3);
    p3=plot([cbf2 cbf2]*aScale,[min(min(min(profBook2))),max(max(max(profBook2)))],'color',[.5 .5 .5],'linestyle','--','linewidth',1);
    p4=plot([.4 .4],[min(min(min(profBook2))),2.2],'color',[1,.5,.5],'linestyle','--','linewidth',1);
    p5=plot([1.6 1.6],[min(min(min(profBook2))),max(max(max(profBook2)))],'color',[1,.5,.5],'linestyle','--','linewidth',1);
    set(gca,'Color',bcolor)
    text(1.05,-.8,'\leftarrowAvg. Cell Boundary','fontsize',16,'color',fcolor)
    
    set(gca,'fontsize',28,'XColor',fcolor,'YColor',fcolor,'YMinorTick','on','LineWidth',2)
    xt = 'Location on Trace (AU)';% input('enter the xaxis label','s');
    yt = 'Displacement (\mum)'; %input('enter the yaxis label','s');
    tt = 'Line-Profile Displacements';%input('enter the title','s');
    le = 'Shear'; %input('enter the legend','s');
    le2 = 'Normal';
    le3 = 'Border';
    xl = xlabel(xt);
    yl = ylabel(yt);
    %tl = title(tt);
    ytickformat('%.1f')
    set(xl, 'fontweight','bold','fontsize',28,'color',fcolor);
    set(yl,'fontweight','bold','fontsize',28,'color',fcolor);
    leg = legend([p1 p2],['\color{' fcolor '}' le],['\color{' fcolor '}' le2],'location','northwest');
    leg.FontSize = 20;
    legend boxoff
    %set(tl,'fontweight','bold','fontsize',title_font_size)
    axis([0 2 min(min(min(profBook2))) max(max(max(profBook2)))])
    
    
    title = strcat(prefix,['All_Data ' fcolor ' on ' bcolor]);
    savefile = [ListPath title];
    export_fig(ProfileOverlays,savefile,'-native');
    
    
    % Normalized Graph
    profBook3 = profBook2;
    for i = 1:size(profBook2,3)
        top(i) = abs(max(max(max(profBook2(2:3,:,i)))));
        profBook3(2:3,:,i) = profBook3(2:3,:,i)/top(i);
    end
    ProfileOverlaysNorm = figure;
    hold on
    for i = 1:size(profBook3,3)
        plot(profBook3(1,:,sortIdx(i))*aScale,profBook3(2,:,sortIdx(i)),'color',[map(i,1:3)])
        plot(profBook3(1,:,sortIdx(i))*aScale,profBook3(3,:,sortIdx(i)),'color',[map(i,1:3)],'linestyle','-.')
    end
    p1 = plot(profBook3(1,:,sortIdx(i))*aScale,mean(profBook3(2,:,sortIdx(:)),3,'omitnan'),'color',fcolor,'linewidth',3);
    p2 = plot(profBook3(1,:,sortIdx(i))*aScale,mean(profBook3(3,:,sortIdx(:)),3,'omitnan'),'color',fcolor,'linestyle','-.','linewidth',3);
    p3 = plot([cbf2 cbf2]*aScale,[min(min(min(profBook3))),1],'color',[.5 .5 .5],'linestyle','--','linewidth',1);
    p4=plot([.4 .4],[min(min(min(profBook2))),.6],'color',[1,.5,.5],'linestyle','--','linewidth',1);
    p5=plot([1.6 1.6],[min(min(min(profBook2))),max(max(max(profBook2)))],'color',[1,.5,.5],'linestyle','--','linewidth',1);
    set(gca,'Color',bcolor)
    text(1.05,-.25,'\leftarrowAvg. Cell Boundary','fontsize',16,'color',fcolor)
    
    set(gca,'fontsize',28,'XColor',fcolor,'YColor',fcolor,'YMinorTick','on','LineWidth',2)
    xt = 'Location on Trace (AU)';% input('enter the xaxis label','s');
    yt = 'Displacement (AU)'; %input('enter the yaxis label','s');
    tt = 'Line-Profile Displacements';%input('enter the title','s');
    le = 'Shear'; %input('enter the legend','s');
    le2 = 'Normal';
    le3 = 'Border';
    xl = xlabel(xt);
    yl = ylabel(yt);
    %tl = title(tt);
    set(xl, 'fontweight','bold','fontsize',28,'color',fcolor);
    set(yl,'fontweight','bold','fontsize',28,'color',fcolor);
    ytickformat('%.2f')
    leg = legend([p1 p2],['\color{' fcolor '}' le],['\color{' fcolor '}' le2],'location','northwest');
    leg.FontSize = 20;
    legend boxoff
    axis([0 2 min(min(min(profBook3))) max(max(max(profBook3)))])
    ytickformat('%.1f')
    %set(tl,'fontweight','bold','fontsize',title_font_size)
    
    
    title = strcat(prefix,['All_Data_Normalized ' fcolor ' on ' bcolor]);
    savefile = [ListPath title];
    export_fig(ProfileOverlaysNorm,savefile,'-native');
end
%%
for i = 1:size(dirList,1)
    cd(dirList{i})
    load('SphereIndent.mat')
    figure
    % Heatmap for XZ
    colorMapRad = colormap(jet(77000));
    colorMapRad(65537:end,:) = [];
    colorMapRad = flipud(colorMapRad);
    MaximumHeatMap = imagesc(radXX(:,1),radYY(:,1),radvq);
    radHeat = MaximumHeatMap.CData;%.*(imageBinary==0);
    
    radHeatNaN = (isnan(radHeat));
    radHeat(isnan(radHeat)) = 0;
    heatScale = (65536/2.4); %(max(max(imageHeat)))
    radHeat = uint16(round(radHeat * heatScale));
    radHeatColor = ind2rgb(radHeat,colorMapRad);
    
    
    scale = (max(max(radYY))-min(min(radYY)))/51;
    addTop = round(min(min(radYY)) / scale);
    addbottom = round((22 - max(max(radYY)))/ scale);
    total = zeros(addTop+addbottom+51,500,3);
    total(addTop:addTop+50,:,:) = radHeatColor(:,:,:);
    
    radHeatColor2 = imresize(total, [100 500]);
    
    radHeatColor3(:,:,:,i) =  radHeatColor2;
    close
    % figure
    % imshow(radHeatColor)
    % hold on
    % text(200,25,[num2str(max(max(radYY))) ' ' num2str(min(min(radYY)))])
    
    figure
    imshow(radHeatColor2)
    hold on
    text(200,25,[num2str(max(max(radYY))) ' ' num2str(min(min(radYY)))])
    
    
end
radHeatColor3(radHeatColor3(:,:,1:3,:)==0) = NaN;
rHC3Mean = mean(radHeatColor3(:,:,:,:),4,'omitnan');
mInd = figure;
imshow(rHC3Mean)
hold on
for i = 1:4
    plot([-1 525],[(500/22)*(i) (500/22)*(i)],'k')
    %text(200,25,[num2str(max(max(radYY))) ' ' num2str(min(min(radYY)))])
end
for i = 1:11
    plot([(2500/50)*(i) (2500/50)*(i)],[-1 105],'k')
    %text(200,25,[num2str(max(max(radYY))) ' ' num2str(min(min(radYY)))])
end
title = ['\Mean Indentation Profile'];
savefile = [ListPath title];
export_fig(mInd,savefile,'-native');
cd(ListPath)

%% Calculate Overlap of Dilated Cell Area and Z-Deformation
%This code is to determine whether 8.125 micron dilation (50pixels@.1625)
%is sufficient to encompass all Z-deformation information
goodData = [1;2;3;5;6;14;17;23;31;34;43;44;47;49;50;48;11;12;19;25;28;36];
clear OverlapPct
for i = 1:size(goodData)
    OverlapPct(i) = nonDefOverlap(dirList{goodData(i)});
end

%% Radial Profiles
clear Radials Maxes
for i = 1:size(dirList,1)
    %try
    [MaxY, MaxX, RadY, RadX, radius] = RadialProfiles(dirList{i,1},prefix);
    radii(i,1) = radius*.1625;
    Maxes(1:size(MaxY,1)) = MaxY;
    Radials(1:size(RadY,1),1:size(RadY,2),i) = RadY;
    load([dirList{i,1} '\Shear Mat Files\ShearXYCutOff.mat'])
    cutOffs(i,1) = xyStd2;
    disp(num2str(i))
    %catch
    %   disp('Radial Profiles Failed')
    %  i
    %end
end



clear RadialsM RadialsM2
for i = 1:size(Radials,2)
    for j = 1:size(Radials,3)
        zeroI = find(Radials(4:end,i,j)<.3,1,'first');
        Radials(zeroI+3:end,i,j) = 0;
    end
end

%%
Radials(Radials == 0) = NaN;
Radials(Radials < mean(cutOffs)/1000) = NaN;
Radials(isnan(Radials)) = .3;
RadialsM = mean(Radials,3,'omitnan');
for i = 1:round(size(RadialsM,2)/2)
    RadialsM2(:,i) = mean([RadialsM(:,i),RadialsM(:,end-(i-1))],2,'omitnan');
end

%Pair up opposing slices (i.e. opposite angle relative to horizon)
for i = 1:round(size(Radials,2)/2)
    pairs(i,1) = i;
    pairs(i,2) = (size(Radials,2)+1)-i;
end

%Create colormap for associating profiles to pairs of slices
colormap = brewermap(size(pairs,1),'Dark2');



clear RadialsStd
for i = 1:size(RadialsM2,2)
    p1(:,:) = Radials(:,pairs(i,1),:);
    p2(:,:) = Radials(:,pairs(i,2),:);
    RadialsStd(:,i) = std([p1,p2],0,2,'omitnan') ;
end

RadialsEr = RadialsStd./sqrt(size(p1,2)*2);

cd(ListPath)
save('RadialData.mat')

%%
MeanProfile = figure;
hold on
Xs  = 0:.1625:(.1625*(size(RadialsM2,1)-1));
for i = 1:size(RadialsM2,2)
    plot(Xs,RadialsM2(:,i),'Color',colormap(i,:))
    %plot(Xs,RadialsM2(:,i)+RadialsStd(:,i),'Color',colormap(i,:))
    %plot(Xs,RadialsM2(:,i)-RadialsStd(:,i),'Color',colormap(i,:))
end
axis([0 25 0 6])
xt = 'Distance From Pattern (\mum)';
yt = 'Magnitude of Displacement (\mum)';
xl = xlabel(xt);
yl = ylabel(yt);
%tl = title(tt);
set(xl, 'fontweight','bold','fontsize',14);
set(yl,'fontweight','bold','fontsize',14);
title = ['\Mean ' prefix ' Plots'];
savefile = [ListPath title];
export_fig(MeanProfile,savefile,'-native');
cd(ListPath)
%%
%CHANGE FONT SIZES HERE
AxisFontSize = 28;
AxisTitleFontSize = 28;
LegendFontSize = 20;

AllProfiles = figure;
set(gcf,'position',[100,100,600,900])
hold on
for j = 1:size(Radials,2)
    [cidx, ~] = find(pairs == j, 1, 'first');
    for i = 1:size(Radials,3)
        plot(Xs,Radials(:,j,i),'Color',colormap(cidx,:))
    end
end
axis([0 25 0 7])
xt = {'Distance From'; 'Pattern (\mum)'};
yt = {'Magnitude of' ;'Displacement (\mum)'};
set(gca,'fontsize',AxisFontSize)
tt = 'Line-Profile Displacements';%input('enter the title','s');
le = 'Shear'; %input('enter the legend','s');
le2 = 'Normal';
le3 = 'Border';
xl = xlabel(xt);
yl = ylabel(yt);
%tl = title(tt);

set(xl, 'fontweight','bold','fontsize',AxisTitleFontSize);
set(yl,'fontweight','bold','fontsize',AxisTitleFontSize);
leg = legend([p1 p2 p3],le,le2,le3,'location','northwest');
leg.FontSize = LegendFontSize;
%set(tl,'fontweight','bold','fontsize',title_font_size)
title = ['\All ' prefix ' Plots'];
savefile = [ListPath title];
export_fig(AllProfiles,savefile,'-native');
cd(ListPath)

%% Plot Averages over individuals
colOptions{1,1} = 'white';
colOptions{2,1} = 'black';
colOptions{1,2} = 'black';
colOptions{2,2} = 'white';

for cI = 1:2 
        fcolor = colOptions{1,cI};
    bcolor = colOptions{2,cI};
    set(0,'defaultfigurecolor',bcolor)
%CHANGE FONT SIZES HERE
AxisFontSize = 40;
AxisTitleFontSize = 40;
LegendFontSize = 30;
Xs  = 0:.1625:(.1625*(size(RadialsM2,1)-1));
AllProfiles2 = figure;
set(gcf,'position',[100,100,1200,650])
hold on
% if strcmp(prefix ,'pill ') == 1
%     rectangle('Position',[0,0,14,7],'FaceColor',[.7 .7 .9020],'EdgeColor',[.7 .7 .9020])
% else
% rectangle('Position',[0,0,8,7],'FaceColor',[.65 .85 1],'EdgeColor',[.9 0 .9])
% end
    
    for j = 1:size(Radials,2)
        
        [cidx, ~] = find(pairs == j, 1, 'first');
        if strcmp(prefix ,'pill ') == 1
            for i = 1:size(Radials,3)
                p1 = plot(Xs,Radials(:,j,i),'Color',colormap(cidx,:),'LineWidth',1.5);
            end
            
        else
            for i = 1:size(Radials,3)
                p1 = plot(Xs,Radials(:,j,i),'Color',[.8 .8 .8],'LineWidth',1.5);
            end
        end
    end
    
    if strcmp(prefix ,'pill ') == 1
        for i = 1:size(RadialsM2,2)
            p2 = plot(Xs,RadialsM2(:,i),'Color',colormap(i,:),'LineWidth',9,'LineStyle','-.');
        end
    else
        
        p2 = plot(Xs,mean(RadialsM2,2,'omitnan')','Color',fcolor,'LineWidth',9,'LineStyle','-.');
              
    end
    
    %plot([0 20],[0.5 0.5],'LineWidth',3,'LineStyle','--','Color','black')
    
    
    
    if strcmp(prefix ,'pill ') == 1
        yticks([0:2:7])
        axis([0 20 0 7])
    else
        yticks([0:3])
        axis([0 20 0 3])
    end
    
    rectangle('Position',[0,0,25,.3],'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7])
    text(9,0.2,'Noise Range','FontSize',20)
    xt = {'Distance From Pattern (\mum)'};
    yt = {'Displacement (\mum)'}; %'Magnitude of' ;
    set(gca,'Color',bcolor)
    set(gca,'fontsize',AxisFontSize)
    set(gca,'fontsize',AxisFontSize,'XColor',fcolor,'YColor',fcolor,'LineWidth',2)
    tt = 'Line-Profile Displacements';%input('enter the title','s');
    le = 'Single Trace'; %input('enter the legend','s');
    le2 = 'Average Trace';
    le3 = 'Border';
    xl = xlabel(xt);
    yl = ylabel(yt);
    %tl = title(tt);
    
    set(xl, 'fontweight','bold','fontsize',AxisTitleFontSize,'color',fcolor);
    set(yl,'fontweight','bold','fontsize',AxisTitleFontSize,'color',fcolor);
    %leg = legend([p1 p2],le,le2,'location','Best');
    %leg.FontSize = LegendFontSize;
    %set(tl,'fontweight','bold','fontsize',title_font_size)
    title = ['\All ' prefix ' Plots 2 '  fcolor ' on ' bcolor];
    savefile = [ListPath title];
    export_fig(AllProfiles2,savefile,'-native');
    cd(ListPath)

end


%%
if strcmp(prefix ,'pill ') == 1
    RadialsM3 = RadialsM2(:,3)';
else
    RadialsM3 = mean(RadialsM2,2,'omitnan')';
end
RadialsM3(isnan(RadialsM3)) = [];
RadialsM4 = RadialsM3-0.5;
Xs2 = [0:0.1625:(0.1625*(size(RadialsM4,2)-1))];
RadM4 = fit(Xs2',RadialsM4','cubicinterp');
figure
plot(RadM4)
if strcmp(prefix ,'pill ') == 1
    fzero(RadM4,14)
else
    fzero(RadM4,7)
end

%%
Xs3 = [0:0.1625:(0.1625*(size(Radials,1)-1))];
if strcmp(prefix ,'pill ') == 1
    for i = 1:size(Radials,3)
        for j = 1:size(Radials,2)
            fun1 = fit(Xs3',Radials(:,j,i)-.5,'cubicinterp');
            zeros(i,j) = fzero(fun1,10);
        end
    end
    
    for i = 1:size(pairs,1)
        RadialsZeros(:,i) = cat(1,zeros(:,pairs(i,1)),zeros(:,pairs(i,2)));
    end
    
    RadialsZS(1,1:3) = mean(RadialsZeros);
    RadialsZS(2,1:3) = std(RadialsZeros);
    save('RadialZerosPill.mat','RadialsZS','zeros')
else
    figure
    hold on
    
    RadialsC = Radials - .5;
    for i = 1:size(RadialsC,3)
        
        
        
        
        
        for j = 1:size(RadialsC,2)
            fun1 = fit(Xs3',RadialsC(:,j,i),'cubicinterp');
            %plot(fun1)
            %xlim([0,10])
            %ylim([-1,10])
            try
                idx = find(RadialsC(:,j,i)>0,1,'last');
                zeros(i,j) = fzero(fun1,[idx]);
                
            catch
                zeros(i,j) = NaN;
            end
            %            plot([zeros(i,j) zeros(i,j)],[0 1])
            
        end
    end
    zeros(zeros>20) = 0;
    zeros(isnan(zeros)) = 0;
    clear RadialsZS
    RadialsZS(1,1) = mean(zeros(:));
    RadialsZS(2,1) = std(zeros(:));
    save('RadialZerosCircle.mat','RadialsZS','zeros')
    
end

%% Max of each group

