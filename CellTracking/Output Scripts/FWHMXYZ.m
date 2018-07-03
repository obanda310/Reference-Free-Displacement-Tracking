%Calculate FWHM dimensions on non-deformed fiducial markers
%This code should provide a quick assessment of marker size using the
%outputs from disp3D.m

%Directory should have 3Ddata.mat file in it.
clear all 
close all
load('3Ddata.mat')

xW = 7;
yW = 7;
zW = 6;

%% X Profiles
%Grab intensity Y as a function of steps in X
figure
hold on
for i = 1:size(m3.ref,1)
    try
clear temp

Xs = (0:raw.dataKey(9,1):(xW*2)*raw.dataKey(9,1));
temp(1,1:2) = round(m3.ref(i,1:2)/raw.dataKey(9,1));
temp(1,3) = round(m3.ref(i,3)/raw.dataKey(10,1));
tempStack = mean(image.RawStack(temp(1,1)-1:temp(1,1)+1,temp(1,2)-xW:temp(1,2)+xW,temp(1,3)),1);
plot(Xs,tempStack)
FWHMs(i,1) = fwhm(Xs,tempStack);
    catch
        i
        FWHMs(i,1) = NaN;
    end
end
FWHMsMean(1,1) = mean(FWHMs(:,1),'omitnan');
FWHMsMean(2,1) = std(FWHMs(:,1),'omitnan');

%% Y Profiles

%Grab intensity Y as a function of steps in X
figure
hold on
for i = 1:size(m3.ref,1)
    try
clear temp

Xs = (0:raw.dataKey(9,1):(xW*2)*raw.dataKey(9,1));
temp(1,1:2) = round(m3.ref(i,1:2)/raw.dataKey(9,1));
temp(1,3) = round(m3.ref(i,3)/raw.dataKey(10,1));
tempStack = mean(image.RawStack(temp(1,1)-yW:temp(1,1)+yW,temp(1,2)-1:temp(1,2)+1,temp(1,3)),2);
plot(Xs,tempStack)
FWHMs(i,2) = fwhm(Xs,tempStack);
     catch
         i
         FWHMs(i,2) = NaN;
     end
end

FWHMsMean(1,2) = mean(FWHMs(:,2),'omitnan');
FWHMsMean(2,2) = std(FWHMs(:,2),'omitnan');



%% Z Profiles
%Grab intensity Y as a function of steps in X
figure
hold on
for i = 1:size(m3.ref,1)
    try
clear temp

Xs = (0:raw.dataKey(10,1):(zW*2)*raw.dataKey(10,1));
temp(1,1:2) = round(m3.ref(i,1:2)/raw.dataKey(9,1));
temp(1,3) = round(m3.ref(i,3)/raw.dataKey(10,1));
tempStack = double(squeeze(mean(mean(image.RawStack(temp(1,1)-1:temp(1,1)+1,temp(1,2)-1:temp(1,2)+1,temp(1,3)-zW:temp(1,3)+zW),2),1)));
plot(Xs,tempStack)
FWHMs(i,3) = fwhm(Xs,tempStack');
     catch
         i
         FWHMs(i,3) = NaN;
     end
end
FWHMs(FWHMs(:,3)<2,3) = NaN;
FWHMsMean(1,3) = mean(FWHMs(:,3),'omitnan');
FWHMsMean(2,3) = std(FWHMs(:,3),'omitnan');

