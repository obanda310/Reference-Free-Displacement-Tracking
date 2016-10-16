function finalMask = createMask(img)
    % Approximate psf and use that psf to deconvolve the image to reduce
    % noise and feature-to-feature blurring. If the psf is known, then
    % normal deconvolution can be performed - must edit code
    
    % Point spread function (kernel) size
    psfSize = 5;
    % Perform blind deconvolution and obtain deconvolved image and
    % approximated psf
    [deconvImg,psf] = deconvblind(img,ones(psfSize));
    % Cross correlate the psf with the deconvolved image to enhance feature
    % brightness (due to high correlation) and lessen inter-feature noise
    % (due to low correlation)
    corrImg = xcorr2(deconvImg,psf);
    
%     % This chunk of code isn't directly useful, but it may give an idea
%     % about the consistency of the spacing of the features in x and y
%     %
%     % Calculate the autocorrelation of the deconvolved image to visualize
%     % the spatial frequencies of the features    
%     acrr = xcorr2(deconvImg);
%     % Determine the indices that correspond to the middle of the
%     % correlation image in x (middle column) and y (middle row)
%     acrrMidRow = ceil(size(acrr,1)/2);
%     acrrMidCol = ceil(size(acrr,2)/2);
%     % Use those indices to consider only the part of the correlation
%     % surface that corresponds to shifts in x and shifts in y (i.e. not
%     % diagonal shifts)
%     acrrX = acrr(acrrMidRow,:);
%     acrrY = acrr(:,acrrMidCol);
%     % Plot the frequency response in x and y to visualize periodicity of
%     % features
%     figure;
%     plot(accrX)
%     hold on
%     plot(accrY)
%     hold off
%     legend('X','Y')

    % Determine size of psf (redundant code, but I think including it makes
    % it easier to read)
    [psfRow,psfCol] = size(psf);
    % The correlation image is larger than the original image because it
    % expands during the correlation process by the size of the psf kernel
    % minus 1. Here, we calculate how many rows and columns we must crop
    % from the correlation image to make it the same size as the original
    % image. We divide by two because we will crop half of the rows that
    % must be cropped from the top edge of the image and the other half 
    % from the bottom edge of the image. Likewise for the columns on the 
    % left and right edges.
    rowShift = (psfRow-1)/2;
    colShift = (psfCol-1)/2;
    corrImgCrop = corrImg(rowShift+1:end-rowShift,colShift+1:end-colShift);
    % Find the peaks in the correlation image to isolate each feature well
    % from neighboring features
    corrPeaks = imregionalmax(corrImgCrop);
    % Dilate the peaks so that the mask covers enough of the features so
    % that we can obtain trustworthy subpixel centroid locations
    corrPeaksDilated = imdilate(corrPeaks,strel('disk',3));
    % Multiply the dilated mask by the original image to recapture the
    % actual shape of the features
    filteredImg = corrPeaksDilated.*img;
    % Generate a new mask that accounts for the true geometry of each
    % feature
    mask = imbinarize(filteredImg,'adaptive');
    % Filter out features that have areas of 5 pixels or less and features
    % that have areas of 30 pixels or more
    finalMask = bwareafilt(mask,[5,30]);
end