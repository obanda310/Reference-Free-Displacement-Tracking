function SphericalIndent(directory)
if nargin ==1
    cd(directory);
end
clear all
close all
%%
%CHANGE FONT SIZES HERE
AxisFontSize = 24;
AxisTitleFontSize = 24;
LegendFontSize = 14;

colOptions{1,1} = 'white';
colOptions{2,1} = 'black';
colOptions{1,2} = 'black';
colOptions{2,2} = 'white';
%% Making the spherical indentation data for FE fitting

load('3Ddata.mat')
tic
for i = 1:size(planesGroups,1)
    pgIdx = find(planesGroups(i,:));
    pgLocs(i,1) = mean(planesLoc2(planesGroups(i,pgIdx)));
end
topPlane = find(pgLocs==min(pgLocs));


    IndPlanes = [1:size(vq3,3)];    

vqtest = vq3(:,:,IndPlanes);
h = fspecial('average', [20 20]);
vqtest(isnan(vqtest)) = 0;
for i = 1:size(vqtest,3)
    vqmask(:,:,i) = bwareaopen(vqtest(:,:,i)<0,50000);
end
vqtestfilt1 = vqtest.*double(vqmask);
for i = 1:size(vqtest,3)
    vqtestfilt2(:,:,i) = filter2(h, vqtestfilt1(:,:,i));
end
vqtestfilt3 = vqtestfilt2*-1;
ShowStack(vqtestfilt3)

%%
figure



imshow(sum(vqN,3,'omitnan'),[])


% for i =1:30
% vqtestproj = sum(vq3,3,'omitnan')*-1;
% vqCO = prctile(vqtestproj(:),i+69);
% vqtestproj(vqtestproj<vqCO) = 0;
% 
% cmap = brewermap(25,'reds');
% 
% %inputs:
% %originalmatrix: the original matrix
% %binarisedmatrix = originalmatrix > threshold; %the thresholded matrix, a logical array
% [rows, cols] = ndgrid(1:size(vqtestproj, 1), 1:size(vqtestproj, 2));
% rc(i) = sum(rows(vqtestproj>0) .* vqtestproj(vqtestproj>0)) / sum(vqtestproj(vqtestproj>0));
% cc(i) = sum(cols(vqtestproj>0) .* vqtestproj(vqtestproj>0)) / sum(vqtestproj(vqtestproj>0));
% 
% 
% hold on
% scatter(cc(i),rc(i),50,'MarkerEdgeColor','r','linewidth',2)
% end


for i =1:5
vqtestproj = sum(vq3,3,'omitnan')*-1;
vqCO = prctile(vqtestproj(:),i+91);
vqtestproj(vqtestproj<vqCO) = 0;

cmap = brewermap(25,'blues');

%inputs:
%originalmatrix: the original matrix
%binarisedmatrix = originalmatrix > threshold; %the thresholded matrix, a logical array
[rows, cols] = ndgrid(1:size(vqtestproj, 1), 1:size(vqtestproj, 2));
rc(i) = sum(rows(vqtestproj>0) .* vqtestproj(vqtestproj>0)) / sum(vqtestproj(vqtestproj>0));
cc(i) = sum(cols(vqtestproj>0) .* vqtestproj(vqtestproj>0)) / sum(vqtestproj(vqtestproj>0));


hold on
scatter(cc(i),rc(i),50,'MarkerEdgeColor',cmap(i,:),'linewidth',2)
end




rcf = round(mean(rc));
ccf = round(mean(cc));
%% Remove planes where center data doesn't exist
clear planeCounts planeCountsCO
%Identify which planes have fewer members than others
for i = 1:size(planesGroups,1)
    planeCounts(i,1) = (nnz(plane.final(:,planesGroups(i,1:nnz(planesGroups(i,:))))));
end
pCMax = max(planeCounts);
planeCounts(:,2) = planeCounts(:,1)/pCMax;

planeCountsCO = planeCounts(:,2)<.8;

pgLocs(planeCountsCO,:) = [];
planesGroups(planeCountsCO,:) = [];
vqtestfilt3(:,:,planeCountsCO) = [];
IndPlanes = [1:size(vqtestfilt3,3)];

%% Create map of distances to center point

distvq = zeros(size(vqN,1),size(vqN,2));
distvq(rcf,ccf)=1;
distvq=round(bwdist(distvq));
figure
imshow(distvq,[])
hold on
scatter(ccf,rcf,50,'MarkerEdgeColor','r','linewidth',2)

%% Create Z-Stack infographic for radial Averaging
figure
imshow(distvq == 500,[])

%%
td = max(distvq(:));

tic
for j = 1:size(pgLocs,1)
    distvqmask = (abs(vqtestfilt3(:,:,j))>0).*distvq;
    for i = 1:td
        distvqA = sum(sum(distvqmask==i));
        radprofile(i,j) = sum(sum((distvq==i).*vqtestfilt3(:,:,j)))/distvqA;
    end
    toc
end

%%
    fcolor = colOptions{2,1};
    bcolor = colOptions{2,2};
set(0,'defaultfigurecolor',bcolor)
RadialProfilesPerZ=figure;
hold on

[~,order] = sort(pgLocs,1);
    
for i = 1:size(order)
    plot([1:1:size(radprofile,1)]*.1625,radprofile(:,order(i)),'DisplayName',num2str(round(pgLocs(order(i),1),1)))
    
end
 
    
    set(gca,'Color',bcolor,'LineWidth',2)
    %Axes, Text, Legends
    ylim([0 ceil(max(radprofile(:)))+1])
    set(gca,'fontsize',AxisFontSize,'XColor',fcolor,'YColor',fcolor,'YMinorTick','on')
    ytickformat('%.1f')
    xt = 'Radial Distance (\mum)';% input('enter the xaxis label','s');
    yt = 'Displacement (\mum)'; %input('enter the yaxis label','s');
    tt = 'Line-Profile Displacements';%input('enter the title','s');
    xl = xlabel(xt);
    yl = ylabel(yt);
    %tl = title(tt);
    
    set(xl, 'fontweight','bold','fontsize',28,'color',fcolor);
    set(yl,'fontweight','bold','fontsize',28,'color',fcolor);
    legend show
    legend boxoff
    leg.FontSize = LegendFontSize;
    %set(tl,'fontweight','bold','fontsize',title_font_size)
    
    
    %Export Image
    title = ['\RadialProfilesPerZ ' fcolor ' on ' bcolor];
    savefile = [filePath title];
    export_fig(RadialProfilesPerZ,savefile,'-native');

%%
clear radprofile2 radprofile3 radprofileXs radprofileZs
for j = 1:size(planesGroups,1)
    planesLocs3(j) = mean(planesLoc2(1,planesGroups(j,1:nnz(planesGroups(j,:)))));
end
planesLocs4 = (planesLocs3(1,IndPlanes));
for i = 1:size(vqtestfilt3,3)
    radprofile2(:,i) = radprofile(:,i) + planesLocs4(i);
end
radprofile3 = double(radprofile(:,1:size(IndPlanes,2)));
radprofileXs = double((1:td)' * raw.dataKey(9,1) *ones(1,size(IndPlanes,2)));
radprofileZs = double(ones(td,1) * planesLocs4  );

[radXX,radYY] = meshgrid(([1:1:td]' * raw.dataKey(9,1)), min(planesLocs4):(max(planesLocs4)-min(planesLocs4))/50:max(planesLocs4));
radXX = double(radXX);
radYY = double(radYY);
radvq = griddata(radprofileXs,radprofileZs,radprofile3,radXX,radYY);
%% Heatmap for XZ
figure
colorMapRad = brewermap(65536,'spectral');
MaximumHeatMap = imagesc(radXX(:,1),radYY(:,1),radvq);
radHeat = MaximumHeatMap.CData;%.*(imageBinary==0);

radHeatNaN = (isnan(radHeat));
radHeat(isnan(radHeat)) = 0;
heatScale = (65536/3); %(max(max(imageHeat)))
radHeat = uint16(round(radHeat * heatScale));
radHeatColor = ind2rgb(radHeat,colorMapRad);
%%
figure
imshow(radHeatColor)
hold on
toc
%%
save('SphereIndent.mat','radHeatColor','radprofile3','radprofileXs','radprofileZs','radYY','radXX','radvq','IndPlanes','planesLocs4','cc','rc','vqtestproj')
%plot([5 10 15 20 25 30 35 40 45 50 55 60 65 70]/xyScale,[1 50 1 50 1 50 1 50 1 50 1 50 1 50] )