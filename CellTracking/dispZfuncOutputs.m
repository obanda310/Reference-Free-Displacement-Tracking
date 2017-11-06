function dispZfuncOutputs(directory)
if nargin ==1
cd(directory);
end
clear all
close all
load('3DnormalData.mat')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE COLOR MAPS OF DISPLACEMENTS IN Z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colorMap = brewermap(65536,'*PuBu');
colorBar1 = zeros(500,25);
range = uint16(round(linspace(65536,1,500)'));
for i = 1:25
    colorBar1(1:500,i) = range;
end
colorBar2 = ind2rgb(colorBar1,colorMap);
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
savefile = [filePath '\HeatMaps\ColorBarZv2.tif'];
export_fig(colorBarSave,savefile,'-native');


cutoff = 0.2; %microns
imageBinaryCombined = ((imageBinaryDilated==0));%+(imageBinary2==0))==0;

%Determine which 'planes' in planesFinal should be the same plane
clear planesGroups
for i = 1:size(planesLoc2,2)
    clear differences
    differences = planesLoc2 - planesLoc2(1,i);
    planesGroups(i,1:size(find(abs(differences)<2),2)) = find(abs(differences)<2)';
end
planesGroups = unique(planesGroups,'rows');
[HeatMapN,vqN] = heatmapZ(r,rDisp3,planesFinal,planesGroups,imageBinaryCombined,xyScale,0,colorMap,0);
[HeatMap,vq1] = heatmapZ(r,rDisp1Filt,planesFinal,planesGroups,imageBinaryCombined,xyScale,noiseCutoff1,colorMap,1);
[HeatMap2,vq2] = heatmapZ(r,rDisp2Filt,planesFinal,planesGroups,imageBinaryCombined,xyScale,noiseCutoff2,colorMap,2);
[HeatMap3,vq3] = heatmapZ(r,rDisp3Filt,planesFinal,planesGroups,imageBinaryCombined,xyScale,noiseCutoff3,colorMap,3);


% Print Data to txt file
% Calculate Useful parameters
[vq1pos,vq1neg,rDisp1PosTotal, rDisp1PosMax, rDisp1NegTotal, rDisp1NegMax]  = vqStats(vq1,planesGroups);
[vq2pos,vq2neg,rDisp2PosTotal, rDisp2PosMax, rDisp2NegTotal, rDisp2NegMax]  = vqStats(vq2,planesGroups);
[vq3pos,vq3neg,rDisp3PosTotal, rDisp3PosMax, rDisp3NegTotal, rDisp3NegMax]  = vqStats(vq3,planesGroups);

% Write to file
cd HeatMaps
planesLocTxt = fopen('Average Z location of planes.txt','wt');
p1Format = 'Plane no. %1.0f is at %.2f microns from the surface \n';
p2Format = 'Plane no. %1.0f has an average fit deviation of %.8f microns from Method 1 fit \n';
p3Format = 'Plane no. %1.0f has an average fit deviation of %.8f microns from Method 2 fit \n';
p4Format = 'Plane no. %1.0f has an average fit deviation of %.8f microns from Method 3 fit \n';
for i = 1:size(planesGroups,1)
    fprintf(planesLocTxt,p1Format,i,mean(planesLoc2(1,planesGroups(i,(planesGroups(i,:)>0)))));
    fprintf(planesLocTxt,p2Format,i,mean(rDisp1MeanPlanes(1,planesGroups(i,(planesGroups(i,:)>0)))));
    fprintf(planesLocTxt,p3Format,i,mean(rDisp2MeanPlanes(1,planesGroups(i,(planesGroups(i,:)>0)))));
    fprintf(planesLocTxt,p4Format,i,mean(rDisp3MeanPlanes(1,planesGroups(i,(planesGroups(i,:)>0)))));
end
fclose(planesLocTxt);
cd(filePath)

%Comma Delimiter Version

for j = 1:3
cd HeatMaps
planesLocComma = fopen(strcat('Planes Data ',num2str(j),'.txt'),'wt');
p1Format = 'Plane\tDistance\tMean\tStd\tPositive Total\tPositive Max\tNegative Total\tNegative Max \n';
p2Format = '%1.0f\t%.2f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f \n';
fprintf(planesLocTxt,p1Format);
for i = 1:size(planesGroups,1)
    fprintf(planesLocTxt,p2Format,i,mean(planesLoc2(1,planesGroups(i,(planesGroups(i,:)>0)))),mean(rDisp2MeanPlanes(1,planesGroups(i,(planesGroups(i,:)>0)))),mean(rDisp2StdPlanes(1,planesGroups(i,(planesGroups(i,:)>0)))),rDisp2PosTotal(i,1),rDisp2PosMax(i,1),rDisp2NegTotal(i,1),rDisp2NegMax(i,1));
end
fclose(planesLocComma);
cd(filePath)
end



%
save 3DnormalData
disp('Scipt has Completed')
% Save data for profile views
%%
filePath = cd;
folderName = 'Profile Data';
mkdir(filePath,folderName)
save('Profile Data\vqZ.mat','vqN','vq1','vq2','vq3','imageTrans','HeatMap','HeatMap2','HeatMap3')