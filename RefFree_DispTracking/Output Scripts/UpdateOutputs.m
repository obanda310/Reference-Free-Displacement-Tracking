function UpdateOutputs(directory)
if nargin == 1
cd(directory)
end

load('3Ddata.mat')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE COLOR MAPS OF DISPLACEMENTS IN Z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    rmdir HeatMaps\3D s
catch
end

mkdir('HeatMaps\3D')
mkdir('HeatMaps\3D\ColorBar')
colorMapZ = single(brewermap(65536,'*PuBu'));
colorBar1 = single(zeros(500,25));
range = uint16(round(linspace(65536,1,500)'));
for i = 1:25
    colorBar1(1:500,i) = range;
end
colorBar2 = ind2rgb(colorBar1,colorMapZ);
for i = 1:10
    colorBar2((i*50)-3:(i*50),13:25,:) = 0;
end
colorBar2(1:3,13:25,:) = 0;

%Save Color Bar Image
close all
colorBarSave = figure;
hold on
imshow(colorBar2);
filePath=cd;
savefile = [filePath '\HeatMaps\3D\ColorBar\ColorBarZ.tif'];
export_fig(colorBarSave,savefile,'-native');

colorMapXY = single(brewermap(65536,'*Spectral'));
colorBar1 = single(zeros(500,25));
range = uint16(round(linspace(65536,1,500)'));
for i = 1:25
    colorBar1(1:500,i) = range;
end
colorBar2 = single(ind2rgb(colorBar1,colorMapXY));
for i = 1:10
    colorBar2((i*50)-3:(i*50),13:25,:) = 0;
end
colorBar2(1:3,13:25,:) = 0;

%Save Color Bar Image
close all
colorBarSave = figure;
hold on
imshow(colorBar2);
filePath=cd;
savefile = [filePath '\HeatMaps\3D\ColorBar\ColorBarShear.tif'];
export_fig(colorBarSave,savefile,'-native');

SE = strel('disk',round(10/.1625));
imageBinaryCombined = image.imgNBds; %imdilate(image.ADil==0,SE);%+(imageBinary2==0))==0;
%

[HeatMapZ,vqZ0,vqXY0] = heatmapZ(r,m3,plane,image,raw.dataKey(9,1),0,colorMapZ,colorMapXY,0);
[HeatMapZ,vqZ,vqXY] = heatmapZ(r,m3,plane,image,raw.dataKey(9,1),1,colorMapZ,colorMapXY,3);


% Print Data to txt file
% Calculate Useful parameters
[vq3pos,vq3neg,m3.PosTotal, m3.PosMax, m3.NegTotal, m3.NegMax]  = vqStats(vqZ,plane.groups);

% Write to file
cd HeatMaps
planesLocTxt = fopen('Average Z location of planes.txt','wt');
p1Format = 'Plane no. %1.0f is at %.2f microns from the surface \n';
p4Format = 'Plane no. %1.0f has an average fit deviation of %.8f microns from Method 3 fit \n';
for i = 1:size(plane.groups,1)
    fprintf(planesLocTxt,p1Format,i,plane.gloc(i));
    fprintf(planesLocTxt,p4Format,i,mean(m3.MeanPlanes(1,plane.groups(i,(plane.groups(i,:)>0)))));
end
fclose(planesLocTxt);
cd(filePath)

%Comma Delimiter Version
cd HeatMaps
printStats(m3,plane.groups,planesLocTxt,plane.loc,'3');
cd(filePath)


%%
disp('Saving Data')
save 3Ddata
%%
% Save data for profile views
filePath = cd;
folderName = 'Profile Data';
mkdir(filePath,folderName)
save('Profile Data\vqZ.mat','vqZ0','vqZ','image','HeatMapN','HeatMap3')

disp(['Script has Completed in ' num2str(toc) ' seconds'])