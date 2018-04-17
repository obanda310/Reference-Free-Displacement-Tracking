function dispEmpty(directory)
if nargin ==1
    cd(directory);
end
%%
clear all
close all
%%
tic
load('3Ddata.mat')

cols2 = r.col;
cols2(:,3) = 1:size(cols2,1);
cols2 = sortrows(cols2,[1 2]);

%%
for i = 1:size(cols2,1)-1
    if cols2(i,1) == cols2(i+1,1) && cols2(i,2)~=0
        cols2(i,4) = cols2(i+1,2)- cols2(i,2);
    else
        cols2(i,4) = NaN;
    end
end
coldist = mean(cols2(:,4),'omitnan');
%%
progressbar('Finding Minimum Intensity Detections')
count = 1;
for i = 1:max(cols2(:,1))
    progressbar(i/max(cols2(:,1)))
    tempcol2 = cols2((cols2(:,1)==i),1:3);
    if size(tempcol2,1)>1
        tempcol2 = sortrows(tempcol2,2);
        for j=1:size(tempcol2,1)-1
            if abs((tempcol2(j+1,2)-tempcol2(j,2))-coldist)<coldist/2
                fun1 = fit([tempcol2(j,2):tempcol2(j+1,2)]',shear.Int(tempcol2(j,2):tempcol2(j+1,2),i),'cubicinterp');
                emptyZ(count,3) = fminsearch(fun1,(tempcol2(j,2)+tempcol2(j+1,2))/2);
                emptyZ(count,4) = i;
                if emptyZ(count,3)>0
                    fun2 = fit([tempcol2(j,2):tempcol2(j+1,2)]',shear.rawX(tempcol2(j,2):tempcol2(j+1,2),i),'cubicinterp');
                    emptyZ(count,1) = feval(fun2,emptyZ(count,3));
                    fun3 = fit([tempcol2(j,2):tempcol2(j+1,2)]',shear.rawY(tempcol2(j,2):tempcol2(j+1,2),i),'cubicinterp');
                    emptyZ(count,2) = feval(fun3,emptyZ(count,3));
                else
                    emptyZ(i,1:3) = NaN;
                end
                count = count+1;                
            end
        end
    end
    
end

emptyZ(:,3) = emptyZ(:,3)*raw.dataKey(10,1);

%%
figure
scatter3(emptyZ(:,1), emptyZ(:,2), emptyZ(:,3),'b')
xlim([0 max(emptyZ(:,1))+5])
ylim([0 max(emptyZ(:,2))+5])
zlim([0 size(shear.rawX,1)*raw.dataKey(10,1)])
hold on
scatter3(r.r(:,1), r.r(:,2), r.r(:,3),'r')

%%
figure
plot([12:24]',shear.Int(12:24,60))
hold on
fun1 = fit([12:24]',shear.Int(12:24,60),'cubicinterp');
plot(fun1)
scatter(fminsearch(fun1,18),feval(fun1,17.6209))
%% Run disp3D protocol to collect zero positions and measure displacements
%Build Raw 3D Object
erraw = RawData3D(res,raw,emptyZ(:,1:3));
erraw = TranscribeR(erraw);

%Build Planes
planeer = PlanesData(erraw);
planeer = nborsPlanesF(planeer,erraw,radXY,radZ);
planeer = growPlanes(planeer,erraw);
[planeer,er] = cleanPlanes(planeer,erraw);
er = RawData3D(res,raw,er);
er = TranscribeR(er);


%Identify Non-Deformed Detections
er = regionCheck(er,image.ADil,raw);

%Build Rows
[er,errows] = buildRows2(er,rowV,planeer.final);

%Format Rows
[errows,errowPlanes,errowPlanesIdx,errowsNDCU,er]=formatRows(errows,planeer,er);

%Calculate plane distance from surface
[erplanesLoc2,erplanesLocFiltList] = placePlanes(er,planeer,Surface2);

%Match to 2D Pillars
[er] = match2D(er,raw,shear,errowsNDCU);

%1st Fit Method
em1 = DispData3D;
em1 = method1fit(em1,er,errows,rowV);
em1 = calcDisp(em1,er,errows,errowsNDCU); %based on fits of non-deformed markers
em1 = calcDispSC(em1,er,shear); %includes shift correction in reference approximation
em1.disp = er.r(:,1:3)-em1.refSC(:,1:3);
em1 = dispStats(em1,planeer,errowPlanesIdx,er,errowsNDCU,erplanesLocFiltList);
em1 = dispNoise(em1,er,erplanesLocFiltList);

%2nd Fit Method
em2 = DispData3D;
em2 = method2fit(em2,em1,er,raw,image.ADil,errows,errowPlanesIdx);
em2 = calcDisp(em2,er,errows,errowsNDCU); %based on fits of non-deformed markers
em2 = calcDispSC(em2,er,shear); %includes shift correction in reference approximation
em2.disp = er.r(:,1:3)-em2.refSC(:,1:3);
em2 = dispStats(em2,planeer,errowPlanesIdx,er,errowsNDCU,erplanesLocFiltList);
em2 = dispNoise(em2,er,erplanesLocFiltList);

%3rd Fit Method
em3 = DispData3D;
em3 = method3fit(em3,em1,em2,errows);
em3 = calcDisp(em3,er,errows,errowsNDCU); %based on fits of non-deformed markers
em3 = calcDispSC(em3,er,shear); %includes shift correction in reference approximation
em3.disp = er.r(:,1:3)-em3.refSC(:,1:3);
em3 = dispStats(em3,planeer,errowPlanesIdx,er,errowsNDCU,erplanesLocFiltList);
em3 = dispNoise(em3,er,erplanesLocFiltList);

% HeatMaps
clear erplanesGroups
for i = 1:size(erplanesLoc2,2)
    clear differences
    differences = erplanesLoc2 - erplanesLoc2(1,i);
    erplanesGroups(i,1:size(find(abs(differences)<2),2)) = find(abs(differences)<2)';
end
erplanesGroups = unique(erplanesGroups,'rows');

[HeatMapN,vqN] = heatmapZ(er.r,em3.disp,planeer.final,erplanesGroups,imageBinaryCombined,image.Borders,raw.dataKey(9,1),0,colorMapZ,colorMapXY,4);
[HeatMap3,vq3] = heatmapZ(er.r,em3.disp,planeer.final,erplanesGroups,imageBinaryCombined,image.Borders,raw.dataKey(9,1),em3.noiseCutoff,colorMapZ,colorMapXY,4);

%%
save('empties.mat','emptyZ','er','em3','planeer')
%%
toc