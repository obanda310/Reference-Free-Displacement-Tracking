% Input variable images is a 3D image stack whose dimensions correspond to
% (rows,columns,z-slice)
function [filtMasks] = getCentroidsStack(images,metadata) %,centroids
    noImgs = size(images,3);
    % Maximum possible intensity value is 2 raised to the "colorDepth" 
    % power - usually, in our case, 16; i.e. we use 16-bit images. Subtract
    % 1 since the range of values begins at 0, not 1.
    maxI = 2^metadata.colorDepth - 1;
    % Normalize image grayscale intensity values so the dynamic range is
    % from 0 to 1 instead of 0 to "maxI"
    images = images/maxI;
    % These next few lines define two standard deviations that will specify
    % the shape of two gaussian smoothing kernels to be used in the
    % for-loop below. These sigma values (i.e. sig1 and sig2) were obtained
    % in code written by Peter Kovesi. These values and the remainder of
    % the code perform a "difference of gaussians" filter operation
    d = 1.2/(metadata.scalingX*1000000);
    sig1 = 1/(1+sqrt(2))*d;
    sig2 = sqrt(2) * sig1;
    % For the last image in the z-stack, find the grayscale intensities 
    % that correspond to the 95th percentile intensity value in each row. 
    % Then, average those values across all rows. This average is equal to 
    % the noiseRng. The last image in the z-stack is chosen because it is
    % black and any non-zero intensity in this image is the result of noise
    noiseRng = mean(prctile(images(:,:,noImgs),95));
    Lap = fspecial('laplacian');
    % Enhance contrast, apply a difference of gaussians filter, then a
    % Laplacian filter to each image in the stack
    lapImages = zeros(size(images));
    for i = 1:noImgs
        thisImg = images(:,:,i);
        % Find highest pixel intensity value in thisImg
        highIn = max(max(max(thisImg)));
        % Adjust intensity values of thisImg such that intensity values
        % less than or equal to noiseRng map to 0 and intensity values
        % greater than or equal to highIn map to 1
        thisImg = imadjust(thisImg,[noiseRng,highIn],[]);
        % Apply a difference of gaussians filter to each image in the
        % stack, with gaussians defined by the sigma values sig1 and sig2
        % determined earlier
        dogImages = 10 * (imgaussfilt(thisImg*maxI, sig2) - ...
                          imgaussfilt(thisImg*maxI, sig1));
        % Apply a Laplacian filter to the 
        lapImages(:,:,i) = imfilter(dogImages,Lap);
    end
    % Determine the 25th percentile pixel intensity value of all images in
    % lapImages, considering only positive pixel intensity values (negative
    % intensity values are introduced after performing the filtering steps)
    % This 25th percentile value will be the threshold used to create a
    % mask for each image such that pixels whose intensities are less than
    % the threshold are set equal to zero, while pixels whose intensities
    % are greater than or equal to the threshold are set equal to one.
    threshold = prctile(lapImages(lapImages>0),25);
    masks = lapImages>threshold;
    % The images in "masks" have some objects that are clearly the result
    % of noise; they are only a few pixels large and not in the expected,
    % ordered grid pattern. This next for-loop eliminates these spurious
    % objects.
    pixelSize = metadata.scalingX*1000000;
    for i = 1:noImgs
        % Use bwareaopen to eliminate objects in "masks" that have an area
        % smaller than 5 pixels
        filtMasks(:,:,i) = bwareaopen(masks(:,:,i),round((0.9/pixelSize)^2));
%         % Use regionprops to find the centroids of all the objects in the
%         % new filtered masks, i.e. filtMasks
%         c = regionprops(filtMasks(:,:,i),'Centroid');
%         % Format the centroid information so that each cell in "centroids"
%         % is a 2-column matrix that has the x and y positions for each
%         % centroid in the image that corresponds to that cell
%         centroids{i} = cat(1,c.Centroid);
    end
end