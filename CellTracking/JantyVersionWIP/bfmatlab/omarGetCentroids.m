% The code I commented out is code that goes unused. If you want to see the
% results of these lines of code anyway, either run this as a script or add
% those variables to the function's output
% Input variable images is a 3D image stack whose dimensions correspond to
% (rows,columns,z-slice)
function [ppImages5,centroids2] = omarGetCentroids(images,metadata)
    noImgs = size(images,3);
    % Maximum possible intensity value is 2 raised to the "colorDepth" 
    % power - usually, in our case, 16; i.e. we use 16-bit images. Subtract
    % 1 since the range of values begins at 0, not 1.
    maxI = 2^metadata.colorDepth - 1;
    % Normalize image grayscale intensity values so the dynamic range is
    % from 0 to 1 instead of 0 to "maxI"
    images = images/maxI;
    
    d = 3;
    sig1 = 1/(1+sqrt(2))*d;
    sig2 = sqrt(2) * sig1;

    % For the last image in the z-stack, find the grayscale intensities 
    % that correspond to the 95th percentile intensity value in each row. 
    % Then, average those values across all rows. This average is equal to 
    % the noiseRng. The last image in the z-stack is chosen because it is
    % black and any non-zero intensity in this image is the result of noise
    noiseRng = mean(prctile(images(:,:,noImgs),95));
%     [mVal,mIdx] = max(prctile(mean(images(:,:,:)),70));
%     refHistImage = images(:,:,mIdx);
%     LoG = fspecial('log',3,.25);
    Lap = fspecial('laplacian');
%     % imadjust(ppImages, [noiseRng;max(max(max(images)))],[0,1];

    parfor i = 1:noImgs
        thisImg = images(:,:,i);
        % Find highest pixel intensity value in thisImg
        highIn = max(max(max(thisImg)));
        % Adjust intensity values of thisImg such that intensity values
        % less than or equal to noiseRng map to 0 and intensity values
        % greater than or equal to highIn map to 1
        thisImg = imadjust(thisImg,[noiseRng,highIn],[]);
        
%         ppImages = imadjust(thisImg/maxI,[noiseRng/maxI,highIn/maxI],[])*maxI;
        ppImages2 = 10 * (imgaussfilt(thisImg*maxI, sig2) - ...
                          imgaussfilt(thisImg*maxI, sig1));
        ppImages3(:,:,i) = imfilter(ppImages2,Lap);
    end
    ppImages4 = (ppImages3>(prctile(ppImages3(ppImages3>0),25)));
    parfor i = 1:noImgs
        ppImages5(:,:,i) = bwareaopen(ppImages4(:,:,i),5);
        c = regionprops(ppImages5(:,:,i),'Centroid');
        centroids2{i} = cat(1,c.Centroid);
    end
end