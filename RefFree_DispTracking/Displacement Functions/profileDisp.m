function [xyq2,zq2,xq,dispProfXY,dispProfZ,heatMap,cell_boundary] = profileDisp(HeatMap,HeatMapXY,vqXY,vqZ,Cent,Max,imageArea,HeatMap3)
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
elseif nargin == 8
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

%%
HeatMapThickProfFig = figure;
imshow(HeatMapThickProf)

%%
savefile = [filePath '\Profile Data\HeatMapThickProfile.tif'];
export_fig(HeatMapThickProfFig,savefile,'-native');

%% Create txt docs for both normalized and original data

Xs3 = Xs/max(Xs);
xq = .01:.01:1;
xyq = interp1(Xs3,Ys(:,1),xq);
zq = interp1(Xs3,Ys(:,2),xq);
dPEFq = interp1(Xs3,dPEF,xq);
%%
figure
plot(xq,xyq);
hold on
plot(xq,zq);
%%
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
dPEFq2 = interp1(xq2,dPEFq,xq);
cell_boundary = find(dPEFq2>0,1,'first');

figure
plot(xq,xyq2);
hold on
plot(xq,zq2);


colors = {'black';'red'};
p = pub_fig([xq',xq'],[xyq2',zq2'],1,'-',colors,'o',1, colors,colors,'ProfThickLinePlotLabeledNormalized.tif');
%%
%This will find the first and last point of the new profile trace found in
%previous section and plot a line and rectangle on a heat map of shear
%deformation to demonstrate where data is being measured at for the
%previous graph

pts4(1,1) = pts2(1,1)+((pts2(1,3)-pts2(1,1))*(newStart/100)); %new x1
pts4(1,2) = pts2(1,2)+((pts2(1,4)-pts2(1,2))*(newStart/100)); %new y1
pts4(1,3) = pts2(1,1)+((pts2(1,3)-pts2(1,1))*(newEnd/100)); %new x2
pts4(1,4) = pts2(1,2)+((pts2(1,4)-pts2(1,2))*(newEnd/100)); %new y2

pts5(1,1) = pts2(1,1)+((pts2(1,3)-pts2(1,1))*(iIdx/100)); %new x1
pts5(1,2) = pts2(1,2)+((pts2(1,4)-pts2(1,2))*(iIdx/100)); %new y1
pts5(1,3) = pts2(1,1)+((pts2(1,3)-pts2(1,1))*(fIdx/100)); %new x2
pts5(1,4) = pts2(1,2)+((pts2(1,4)-pts2(1,2))*(fIdx/100)); %new y2

pts6(1,1) = pts2(1,1)-((pts2(1,3)-pts2(1,1))*100); %new x1
pts6(1,2) = pts2(1,2)-((pts2(1,4)-pts2(1,2))*100); %new y1
pts6(1,3) = pts2(1,1)+((pts2(1,3)-pts2(1,1))*100); %new x2
pts6(1,4) = pts2(1,2)+((pts2(1,4)-pts2(1,2))*100); %new y2

% heatMap2 = double(zeros(size(HeatMap,1),size(HeatMap,2)));
% heatMap3 = insertShape(heatMap2,'Line',pts2);
% heatMap2 = insertShape(heatMap2,'Line',pts4,'color','green');
% dispProfXY= improfile(vqXY,pts(1,:),pts(2,:),'bilinear');
% dispProfZ= improfile(vqZ,pts(1,:),pts(2,:),'bilinear');
% heatMap2 = im2bw(heatMap2,0.01);
% heatMap3 = im2bw(heatMap3,0.01);
% for i = 1:3
%     if i == 1 ||  i == 3
%     heatMapColor2(:,:,i) = heatMap2*-2;
%     heatMapColor3(:,:,i) = heatMap3;
%     else
%     heatMapColor2(:,:,i) = heatMap2;
%     heatMapColor3(:,:,i) = heatMap3;
%     end
% end

% heatMapColor2 = double(heatMapColor2);
% heatMapFinal2 = heatMapColor2  + HeatMap; %+ heatMapColor3

%%
ExpHeat= figure; 
imshow(HeatMap)
hold on
% Plot the numerical values
fSexp = 20;
fSexp2 = 12;
fSexp3 = 16;
nvX = 18;
nvC = [1,1,1];


scatter(pts4(1,1),pts4(1,2),200,'x','MarkerEdgeColor',nvC)
text(pts4(1,1)+nvX,pts4(1,2),'0.0','color',nvC,'fontsize',fSexp2)

scatter(pts5(1,1),pts5(1,2),200,'x','MarkerEdgeColor',nvC)
text(pts5(1,1)+nvX,pts5(1,2),'0.4','color','r','fontsize',fSexp2)

scatter(mean([pts4(1,1),pts5(1,3)]),mean([pts4(1,2),pts5(1,4)]),200,'x','MarkerEdgeColor',nvC)
text(mean([pts4(1,1) pts5(1,3)])+nvX,mean([pts4(1,2),pts5(1,4)]),'0.8','color',nvC,'fontsize',fSexp2)

scatter(mean([pts4(1,3),pts5(1,1)]),mean([pts4(1,4),pts5(1,2)]),200,'x','MarkerEdgeColor',nvC)
text(mean([pts4(1,3) pts5(1,1)])+nvX,mean([pts4(1,4),pts5(1,2)]),'1.2','color',nvC,'fontsize',fSexp2)

scatter(pts5(1,3),pts5(1,4),200,'x','MarkerEdgeColor',nvC)
text(pts5(1,3)+nvX,pts5(1,4),'1.6','color','r','fontsize',fSexp2)

scatter(pts4(1,3),pts4(1,4),200,'x','MarkerEdgeColor',nvC)
text(pts4(1,3)+nvX,pts4(1,4),'2.0','color',nvC,'fontsize',fSexp2)


% Plot Center and Max Traction
scatter(Cent(1,1),Cent(1,2),50,'o','w')
text(Cent(1,1)+175,Cent(1,2)-25,{'Cell','Centroid'},'color','w','fontsize',fSexp,'HorizontalAlignment','center')
plot([Cent(1,1),Cent(1,1)+120],[Cent(1,2),Cent(1,2)-25],'LineWidth',1,'Color','w')

scatter(Max(1,1),Max(1,2),50,'o','w')
text(Max(1,1)+125,Max(1,2),{'Max','Shear'},'color','w','fontsize',fSexp,'HorizontalAlignment','center')
plot([Max(1,1),Max(1,1)+80],[Max(1,2),Max(1,2)],'LineWidth',1,'Color','w')



%Plot the original trace through center and max Traction
plot([pts6(1,1),pts6(1,3)],[pts6(1,2),pts6(1,4)],'LineStyle','--','LineWidth',1,'Color',[.7,.7,.7])%,'g')
%Plot the profile line
plot([pts4(1,1),pts4(1,3)],[pts4(1,2),pts4(1,4)],'LineWidth',3,'Color','w')%,'g')

%Plot 20th Percentile Tractions
scatter(pts5(1,1),pts5(1,2),1000,'.','r')
text(pts5(1,1)-125,pts5(1,2)+50,{'20% Shear',' (Proximal)'},'color','w','HorizontalAlignment','center','fontsize',fSexp3)
plot([pts5(1,1),pts5(1,1)-50],[pts5(1,2),pts5(1,2)+50],'LineWidth',1,'Color','w')

scatter(pts5(1,3),pts5(1,4),1000,'.','r')
text(pts5(1,3)-125,pts5(1,4),{'20% Shear',' (Distal)'},'color','w','HorizontalAlignment','center','fontsize',fSexp3)
plot([pts5(1,3),pts5(1,3)-50],[pts5(1,4),pts5(1,4)],'LineWidth',1,'Color','w')

savefile = [filePath '\Profile Data\HeatMapProfileExplanatory.tif'];
export_fig(ExpHeat,savefile,'-native');
%%
ExpHeat2 = figure; 
imshow(HeatMapXY)
hold on
% Plot the numerical values
fSexp = 20;
fSexp2 = 12;
fSexp3 = 16;
nvX = 18;
nvC = [1,1,1];
dashL = [0,0,0];
twentyPct = 'r';

%The plus and minus vectors to make tick marks.
perp2= [-1*(pts4(1,4)-pts4(1,2)) pts4(1,3)-pts4(1,1)];
perp2Scale = sqrt((40)/(perp2(1,1)^2+perp2(1,2)^2));
perp2Scale2 = sqrt((200)/(perp2(1,1)^2+perp2(1,2)^2));
perp2Unit(1,1) = perp2(1,1)*perp2Scale;
perp2Unit(1,2) = perp2(1,2)*perp2Scale;
perp2Unit(1,3) = perp2(1,1)*perp2Scale2;
perp2Unit(1,4) = perp2(1,2)*perp2Scale2;



%Plot the original trace through center and max Traction
plot([pts6(1,1),pts6(1,3)],[pts6(1,2),pts6(1,4)],'LineStyle','--','LineWidth',1,'Color',dashL)%,'g')

%Plot the profile line
plot([pts4(1,1),pts4(1,3)],[pts4(1,2),pts4(1,4)],'LineWidth',3,'Color','w')%,'g')

%Plot the tick markers on line
% scatter(pts4(1,1),pts4(1,2),200,'x','MarkerEdgeColor',nvC)
% scatter(pts5(1,1),pts5(1,2),200,'x','MarkerEdgeColor',nvC)
% scatter(mean([pts4(1,1),pts5(1,3)]),mean([pts4(1,2),pts5(1,4)]),200,'x','MarkerEdgeColor',nvC)
% scatter(mean([pts4(1,3),pts5(1,1)]),mean([pts4(1,4),pts5(1,2)]),200,'x','MarkerEdgeColor',nvC)
% scatter(pts5(1,3),pts5(1,4),200,'x','MarkerEdgeColor',nvC)
% scatter(pts4(1,3),pts4(1,4),200,'x','MarkerEdgeColor',nvC)

plot([pts4(1,1)+perp2Unit(1,3) pts4(1,1)-perp2Unit(1,3)],[pts4(1,2)+perp2Unit(1,4) pts4(1,2)-perp2Unit(1,4)],'Color',nvC,'LineWidth',2)
plot([pts5(1,1)+perp2Unit(1,1) pts5(1,1)-perp2Unit(1,1)],[pts5(1,2)+perp2Unit(1,2) pts5(1,2)-perp2Unit(1,2)],'Color',twentyPct,'LineWidth',3)
plot([mean([pts4(1,1),pts5(1,3)])+perp2Unit(1,1) mean([pts4(1,1),pts5(1,3)])-perp2Unit(1,1)],[mean([pts4(1,2),pts5(1,4)])+perp2Unit(1,2) mean([pts4(1,2),pts5(1,4)])-perp2Unit(1,2)],'Color',nvC)
plot([mean([pts4(1,3),pts5(1,1)])+perp2Unit(1,1) mean([pts4(1,3),pts5(1,1)])-perp2Unit(1,1)],[mean([pts4(1,4),pts5(1,2)])+perp2Unit(1,2) mean([pts4(1,4),pts5(1,2)])-perp2Unit(1,2)],'Color',nvC)
plot([pts5(1,3)+perp2Unit(1,1) pts5(1,3)-perp2Unit(1,1)],[pts5(1,4)+perp2Unit(1,2) pts5(1,4)-perp2Unit(1,2)],'Color',twentyPct,'LineWidth',3)
plot([pts4(1,3)+perp2Unit(1,3) pts4(1,3)-perp2Unit(1,3)],[pts4(1,4)+perp2Unit(1,4) pts4(1,4)-perp2Unit(1,4)],'Color',nvC,'LineWidth',2)

% Plot Center and Max Traction
scatter(Cent(1,1),Cent(1,2),1000,'.','MarkerEdgeColor',dashL )
scatter(Max(1,1),Max(1,2),1000,'.','MarkerEdgeColor',dashL )

% %Plot 20th Percentile Tractions
% scatter(pts5(1,1),pts5(1,2),250,'x','black','LineWidth',3)
% scatter(pts5(1,3),pts5(1,4),250,'x','black','LineWidth',3)

savefile = [filePath '\Profile Data\HeatMapProfileExplanatory2.tif'];
export_fig(ExpHeat2,savefile,'-native');

%%
NormalProfile = figure; 
imshow(HeatMap3)
hold on
%The plus and minus vectors to make tick marks.
perp2= [-1*(pts4(1,4)-pts4(1,2)) pts4(1,3)-pts4(1,1)];
perp2Scale = sqrt((40)/(perp2(1,1)^2+perp2(1,2)^2));
perp2Scale2 = sqrt((200)/(perp2(1,1)^2+perp2(1,2)^2));
perp2Unit(1,1) = perp2(1,1)*perp2Scale;
perp2Unit(1,2) = perp2(1,2)*perp2Scale;
perp2Unit(1,3) = perp2(1,1)*perp2Scale2;
perp2Unit(1,4) = perp2(1,2)*perp2Scale2;

%Plot the original trace through center and max Traction
plot([pts6(1,1),pts6(1,3)],[pts6(1,2),pts6(1,4)],'LineStyle','--','LineWidth',1,'Color',dashL)%,'g')

%Plot the profile line
plot([pts4(1,1),pts4(1,3)],[pts4(1,2),pts4(1,4)],'LineWidth',3,'Color','w')%,'g')

%Plot the tick markers on line
% scatter(pts4(1,1),pts4(1,2),200,'x','MarkerEdgeColor',nvC)
% scatter(pts5(1,1),pts5(1,2),200,'x','MarkerEdgeColor',nvC)
% scatter(mean([pts4(1,1),pts5(1,3)]),mean([pts4(1,2),pts5(1,4)]),200,'x','MarkerEdgeColor',nvC)
% scatter(mean([pts4(1,3),pts5(1,1)]),mean([pts4(1,4),pts5(1,2)]),200,'x','MarkerEdgeColor',nvC)
% scatter(pts5(1,3),pts5(1,4),200,'x','MarkerEdgeColor',nvC)
% scatter(pts4(1,3),pts4(1,4),200,'x','MarkerEdgeColor',nvC)

plot([pts4(1,1)+perp2Unit(1,3) pts4(1,1)-perp2Unit(1,3)],[pts4(1,2)+perp2Unit(1,4) pts4(1,2)-perp2Unit(1,4)],'Color',nvC,'LineWidth',2)
plot([pts5(1,1)+perp2Unit(1,1) pts5(1,1)-perp2Unit(1,1)],[pts5(1,2)+perp2Unit(1,2) pts5(1,2)-perp2Unit(1,2)],'Color',twentyPct,'LineWidth',3)
plot([mean([pts4(1,1),pts5(1,3)])+perp2Unit(1,1) mean([pts4(1,1),pts5(1,3)])-perp2Unit(1,1)],[mean([pts4(1,2),pts5(1,4)])+perp2Unit(1,2) mean([pts4(1,2),pts5(1,4)])-perp2Unit(1,2)],'Color',nvC)
plot([mean([pts4(1,3),pts5(1,1)])+perp2Unit(1,1) mean([pts4(1,3),pts5(1,1)])-perp2Unit(1,1)],[mean([pts4(1,4),pts5(1,2)])+perp2Unit(1,2) mean([pts4(1,4),pts5(1,2)])-perp2Unit(1,2)],'Color',nvC)
plot([pts5(1,3)+perp2Unit(1,1) pts5(1,3)-perp2Unit(1,1)],[pts5(1,4)+perp2Unit(1,2) pts5(1,4)-perp2Unit(1,2)],'Color',twentyPct,'LineWidth',3)
plot([pts4(1,3)+perp2Unit(1,3) pts4(1,3)-perp2Unit(1,3)],[pts4(1,4)+perp2Unit(1,4) pts4(1,4)-perp2Unit(1,4)],'Color',nvC,'LineWidth',2)

% Plot Center and Max Traction
scatter(Cent(1,1),Cent(1,2),1000,'.','MarkerEdgeColor',dashL )
scatter(Max(1,1),Max(1,2),1000,'.','MarkerEdgeColor',dashL )


% %Plot 20th Percentile Tractions
% scatter(pts5(1,1),pts5(1,2),250,'x','black','LineWidth',3)
% scatter(pts5(1,3),pts5(1,4),250,'x','black','LineWidth',3)
savefile = [filePath '\Profile Data\TransmittedProfileExplanatory2.tif'];
export_fig(NormalProfile,savefile,'-native');

%%
save('profileDispOutput.mat')
end