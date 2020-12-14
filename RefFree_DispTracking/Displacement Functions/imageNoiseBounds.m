%improving the interpolation and noise filtering scheme.
function [imgNBds6] = imageNoiseBounds(r,m3,image,scale)
%%
imgNBds = zeros(size(image.Black,1),size(image.Black,2));

%%

%Using forced thresholds instead (because a reliable statistical method is not coded)
xyCO = .2; %microns
zCO = .2; %microns

dispXY = sqrt(m3.disp(:,1).^2+m3.disp(:,2).^2);

SE = strel('disk',round(2/scale));
SEarea = round(3.141*3*((round(2/scale))^2));
for i = 1:size(r,1)
    try
        if (abs(m3.disp(i,3)) > zCO  || dispXY(i,1) > xyCO) && round(m3.refSC(i,2)/scale)<=size(imgNBds,1) && round(m3.refSC(i,1)/scale)<=size(imgNBds,2)
        imgNBds(round(m3.refSC(i,2)/scale),round(m3.refSC(i,1)/scale)) = 1;
        end
    catch
        %disp('problems!')
        
    end
end
imgNBds2 = imdilate(imgNBds,SE);
%disp(size(imgNBds2))
imgNBds3 = bwareaopen(imgNBds2==1,(SEarea));
%disp(size(imgNBds3))
imgNBds4 = imdilate(bwareaopen(imgNBds3==1,(SEarea*4)),SE);
%disp(size(imgNBds4))
imgNBds5 = bwareaopen((imgNBds4 | imgNBds3)==0,SEarea)==0;
%disp(size(imgNBds5))
imgNBds5 = imclose(imgNBds5,SE);
imgNBds6 = bwareaopen(imgNBds5==1,(SEarea*4));
% figure
% imshow(imgNBds6,[])