clear all
close all
scaleZ = 0.4;
Z1 = getImages;
loc = find(Z1 == max(max(max(Z1))));
[Z2,Z3,Z4] = ind2sub([size(Z1,1),size(Z1,2), size(Z1,3)],loc);
Z5= Z1(round(mean(Z2)),round(mean(Z3)),:);
Z5 = squeeze(Z5);
xiZ=1:1:size(Z5,1);
xqZ=1:.1:size(Z5,1);
scaleFacZ = size(xiZ,2)/size(xqZ,2);
scaleZ2 = scaleZ*scaleFacZ;
vqZ = interp1(xiZ,Z5,xqZ,'PCHIP');
imshow(vqZ,[])
FWHMz = vqZ(vqZ>=(max(vqZ))/2);
FWHMz2 = size(FWHMz,2)*scaleZ2;

