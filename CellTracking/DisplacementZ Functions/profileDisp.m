function [xyq2,zq2,xq,dispProfXY,dispProfZ,heatMap,cell_boundary] = profileDisp(HeatMap,HeatMapXY,vqXY,vqZ,Cent,Max,imageArea)
filePath = cd;
vqXY(isnan(vqXY))=0;
vqZ(isnan(vqZ))=0;
HeatMap = (HeatMap/3 + HeatMapXY/1.5);
imshow(HeatMap,[]);
HeatMapProfLine = figure;
if nargin == 5
pts=readPoints(HeatMap,2);
for i = 1:size(pts,2)
    pts2(1,1+(2*(i-1))) = pts(1,i);
    pts2(1,2+2*(i-1)) = pts(2,i);
end
elseif nargin == 7
    plusX = (Max(1,1) - Cent(1,1))/1.333;
    plusY = (Max(1,2) - Cent(1,2))/1.333;
    pts(1,1) = Cent(1,1)-plusX;
    pts(2,1) = Cent(1,2)-plusY;
    pts(1,2) = Max(1,1)+plusX;
    pts(2,2) = Max(1,2)+plusY;   
    pts2(1,1) = Cent(1,1)-plusX;
    pts2(1,2) = Cent(1,2)-plusY;
    pts2(1,3) = Max(1,1)+plusX;
    pts2(1,4) = Max(1,2)+plusY;
end

heatMap = double(zeros(size(HeatMap,1),size(HeatMap,2)));
heatMap = insertShape(heatMap,'Line',pts2);
dispProfXY= improfile(vqXY,pts(1,:),pts(2,:),'bilinear');
dispProfZ= improfile(vqZ,pts(1,:),pts(2,:),'bilinear');
heatMap = im2bw(heatMap,0.01);
for i = 1:3
    heatMapColor(:,:,i) = heatMap;
end
heatMapColor = double(heatMapColor);
heatMapFinal = heatMapColor + HeatMap;
figure 
imshow(heatMapFinal)


ProfLineBinary = figure;
imshow(heatMap)

ProfLinePlot = figure;
plot(1:1:size(dispProfXY),dispProfXY)
hold on
plot(1:1:size(dispProfZ),dispProfZ)

Xs = (.1625:.1625:size(dispProfXY)*.1625)';
Xs2 = [Xs,Xs];
Ys = [dispProfXY,dispProfZ];
savefile = [filePath '\Profile Data\HeatMapwLine.tif'];
export_fig(HeatMapProfLine,savefile,'-native');
savefile = [filePath '\Profile Data\ProfLineBinary.tif'];
export_fig(ProfLineBinary,savefile,'-native');
savefile = [filePath '\Profile Data\ProfLinePlots.tif'];
export_fig(ProfLinePlot,savefile,'-native');


colors = {'black';'red'};
p = pub_fig(Xs2,Ys,1,'-',colors,'o',1, colors,colors,'ProfLinePlotLabeled.tif');



%% Create Thicker Profiles
%Code should get line profiles for several lines parallel to the first line
%and average them for "potentially" cleaner data
perp= [-1*(pts(2,1)-pts(2,2)) pts(1,1)-pts(1,2)];
perpScale = sqrt((1)/(perp(1,1)^2+perp(1,2)^2));
perpUnit(1,1) = perp(1,1)*perpScale;
perpUnit(1,2) = perp(1,2)*perpScale;
profWidth = 20;
profVector = perpUnit*profWidth;


pts3 = zeros(1,8);
pts3(1,1:2) = pts2(1,1:2)+profVector;
pts3(1,3:4) = pts2(1,1:2)-profVector;
pts3(1,5:6) = pts2(1,3:4)-profVector;
pts3(1,7:8) = pts2(1,3:4)+profVector;

for i = 1:profWidth
    
    %Shear Data
    ptsTemp = pts;
    ptsTemp(1,1:2) = ptsTemp(1,1:2)+perpUnit(1,1)*i;
    ptsTemp(2,1:2) = ptsTemp(2,1:2)+perpUnit(1,2)*i;
    clear profTemp
    profTemp = improfile(vqXY,ptsTemp(1,:),ptsTemp(2,:),'bilinear');
    current = 2*i-1;
    dispProfXYThck(1:size(profTemp,1),current)= profTemp(:,1);
    ptsTemp = pts;
    ptsTemp(1,1:2) = ptsTemp(1,1:2)-perpUnit(1,1)*i;
    ptsTemp(2,1:2) = ptsTemp(2,1:2)-perpUnit(1,2)*i;
    clear profTemp
    profTemp = improfile(vqXY,ptsTemp(1,:),ptsTemp(2,:),'bilinear');
    current = 2*i;
    dispProfXYThck(1:size(profTemp,1),current)= profTemp(:,1);
    
    
    %Normal Data
    ptsTemp = pts;
    ptsTemp(1,1:2) = ptsTemp(1,1:2)+perpUnit(1,1)*i;
    ptsTemp(2,1:2) = ptsTemp(2,1:2)+perpUnit(1,2)*i;
    clear profTemp
    profTemp = improfile(vqZ,ptsTemp(1,:),ptsTemp(2,:),'bilinear');
    current = 2*i-1;
    dispProfZThck(1:size(profTemp,1),current)= profTemp(:,1);
    ptsTemp = pts;
    ptsTemp(1,1:2) = ptsTemp(1,1:2)-perpUnit(1,1)*i;
    ptsTemp(2,1:2) = ptsTemp(2,1:2)-perpUnit(1,2)*i;
    clear profTemp
    profTemp = improfile(vqZ,ptsTemp(1,:),ptsTemp(2,:),'bilinear');
    current = 2*i;
    dispProfZThck(1:size(profTemp,1),current)= profTemp(:,1);
    
    %Cell Boundary
    ptsTemp = pts;
    clear profTemp
    profTemp = improfile(imageArea,ptsTemp(1,:),ptsTemp(2,:),'bilinear');
    dispProfEdge(1:size(profTemp,1),1)= profTemp(:,1);
end
dispProfXYThckFinal = mean(dispProfXYThck,2);
dispProfZThckFinal = mean(dispProfZThck,2);
dPEF = mean(dispProfEdge,2);
Ys = [dispProfXYThckFinal,dispProfZThckFinal];


colors = {'black';'red'};
p = pub_fig(Xs2,Ys,1,'-',colors,'o',1, colors,colors,'ProfThickLinePlotLabeled.tif');

HeatMapThickProf = insertShape(heatMapFinal,'FilledPolygon',pts3);
HeatMapThickProfFig = figure;
imshow(HeatMapThickProf)
savefile = [filePath '\Profile Data\HeatMapThickProfile.tif'];
export_fig(HeatMapThickProfFig,savefile,'-native');

%% Create txt docs for both normalized and original data

Xs3 = Xs/max(Xs);
xq = .01:.01:1;
xyq = interp1(Xs3,Ys(:,1),xq);
zq = interp1(Xs3,Ys(:,2),xq);
dPEFq = interp1(Xs3,dPEF,xq);
figure
plot(xq,xyq);
hold on
plot(xq,zq);

% Shift curve to be centered on max shear magnitude
pkTails = .2;%Threshold for finding shear peak tails

mIdx = find(xyq == max(xyq));
iIdx = find(round((xyq(1:mIdx)/max(xyq))/pkTails) == 0,1,'last');
fIdx = mIdx + find(round((xyq(mIdx:end)/max(xyq))/pkTails) == 0,1,'first');
pkWidth = fIdx-iIdx;
fWidth = pkWidth*3;
iIdx2 = iIdx-pkWidth/3;
fIdx2 = fIdx+pkWidth/3;
iloc = iIdx2*0.01;
floc = fIdx2*0.01;
if floc>1
    addZero = round(floc-1,2)/0.01;
    xq3 = .01:.01:round(floc,2);
    xyq = cat(2,xyq,zeros(1,addZero));
    zq = cat(2,zq,zeros(1,addZero));
    dPEFq = cat(2,dPEFq,zeros(1,addZero));
else
    xq3 = xq;
end
newStart = find(abs(xq3-iloc) == min(abs(xq3-iloc)));
newEnd = find(abs(xq3-floc) == min(abs(xq3-floc)));
xq2=(xq3-(xq3(newStart)-0.01))/(xq3(newEnd)-(xq3(newStart)-0.01));
xyq2 = interp1(xq2,xyq,xq);
zq2 = interp1(xq2,zq,xq);
dPEFq2 = interp1(xq2,dPEFq,xq)
cell_boundary = find(dPEFq2>0,1,'first');


plot(xq,xyq2);
hold on
plot(xq,zq2);


colors = {'black';'red'};
p = pub_fig([xq',xq'],[xyq2',zq2'],1,'-',colors,'o',1, colors,colors,'ProfThickLinePlotLabeledNormalized.tif');

end