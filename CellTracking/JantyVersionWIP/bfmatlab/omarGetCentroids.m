% The code I commented out is code that goes unused. If you want to see the
% results of these lines of code anyway, either run this as a script or add
% those variables to the function's output
% Input variable images is a 3D image stack whose dimensions correspond to
% (rows,columns,z-slice)
function [ppImages5,centroids2] = omarGetCentroids(images,metadata)
    noImgs = size(images,3);

    d = 3;
    sig1 = 1/(1+sqrt(2))*d;
    sig2 = sqrt(2) * sig1;

    noiseRng = mean(prctile(images(:,:,noImgs),95));
%     [mVal,mIdx] = max(prctile(mean(images(:,:,:)),70));
%     refHistImage = images(:,:,mIdx);
%     LoG = fspecial('log',3,.25);
    Lap = fspecial('laplacian');
%     % imadjust(ppImages, [noiseRng;max(max(max(images)))],[0,1];

    % Maximum possible intensity value is 2 raised to the "colorDepth" 
    % power - usually, in our case, 16; i.e. we use 16-bit images. Subtract
    % 1 since the range of values begins at 0, not 1.
    maxI = 2^metadata.colorDepth - 1;
    parfor i = 1:noImgs
        thisImg = images(:,:,i);
        highIn = max(max(max(thisImg)));
%         ppImages = imadjust(thisImg/maxI,[noiseRng/maxI,highIn/maxI],[])*maxI;
        ppImages2 = ( ...
            imgaussfilt( ...
                imadjust( ...
                    thisImg/maxI, ...
                    [noiseRng/maxI,highIn/maxI], ...
                    []) * ...
                maxI, ...
                sig2) - ...
            imgaussfilt( ...
                imadjust( ...
                    thisImg/maxI, ...
                    [noiseRng/maxI,highIn/maxI], ...
                    []) * ...
                maxI, ...
                sig1) ...
            )*10;
        ppImages3(:,:,i) = imfilter(ppImages2,Lap);
    end
    ppImages4 = (ppImages3>(prctile(ppImages3(ppImages3>0),25)));
    parfor i = 1:noImgs
        ppImages5(:,:,i) = bwareaopen(ppImages4(:,:,i),5);
        c = regionprops(ppImages5(:,:,i),'Centroid');
        centroids2{i} = cat(1,c.Centroid);
    end
end