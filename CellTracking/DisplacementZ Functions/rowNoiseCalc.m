function [noiseMean1,noiseStd1,noiseCutoff1,rDisp1Filt,rDisp1PF] = rowNoiseCalc(r,rDisp1,planesLocFiltList,rNDC)
%% Filter out noise in Displacements Method 2
noiseMean1 = mean(abs(rDisp1(intersect(planesLocFiltList(:,1),rNDC),3)));
noiseStd1 = std((rDisp1(intersect(planesLocFiltList(:,1),rNDC),3)));
noiseCutoff1 = 2*noiseStd1;


% Use noiseCutoff to filter data
rDisp1Filt = rDisp1;
rDisp1Filt(abs(rDisp1Filt)<noiseCutoff1) = 0;
rDisp1PF = find(abs(rDisp1(:,3))<noiseCutoff1+noiseStd1 & abs(rDisp1(:,3))>noiseCutoff1);
clear rNbor
radXY = 2.5;
radZ = .75;
for i = 1:size(r,1)
    topX = r(i,1)+ radXY;
    botX = r(i,1)- radXY;
    topY = r(i,2)+ radXY;
    botY = r(i,2)- radXY;
    topZ = r(i,3)+ radZ;
    botZ = r(i,3)- radZ;
    rNbor(i,1:size(find(r(:,1)<topX & r(:,1)>botX & r(:,2)<topY & r(:,2)>botY& r(:,3)<topZ & r(:,3)>botZ))) = find(r(:,1)<topX & r(:,1)>botX & r(:,2)<topY & r(:,2)>botY& r(:,3)<topZ & r(:,3)>botZ);
end

for i = 1:size(rDisp1PF,1)
    A = rDisp1PF(i,1);
    rDisp1PF(i,2)= rDisp1(rDisp1PF(i,1),3);
    rDisp1PF(i,3) = mean(rDisp1Filt(rNbor(A,rNbor(A,:)~=0&(rNbor(A,:)~=A)),3));
    rDisp1PF(i,4) = std(rDisp1Filt(rNbor(A,rNbor(A,:)~=0&(rNbor(A,:)~=A)),3));
    %filt 2 std from mean
    if pdist([rDisp1PF(i,2);rDisp1PF(i,3)]) > ((rDisp1PF(i,4)))
        rDisp1PF(i,5) = 0;
    else
        rDisp1PF(i,5) = 1;
    end

    if rDisp1PF(i,5) == 0
        rDisp1Filt(rDisp1PF(i,1),:) = 0;
    end
end
end