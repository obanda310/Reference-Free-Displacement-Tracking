% Input mask should be the result of imregional max or any other mask whose
% features are one pixel in size
function centroids = getCentroids2d(mask)
    % Convert the coordinates of mask features (i.e. where the mask image
    % equals one, not zero) to subscript values rowData and colData for use
    % in subpix2d Peter Kovesi function
    [rowData,colData] = ind2sub(size(mask),find(mask==1));
    % Use mask centroids as intial guesses for subpix2d function to find
    % sub-pixel coordinates for the centroids
    [centroidRow,centroidCol] = subpix2d(rowData,colData,mask);
    % Output centroidRow and centroidCol information as an N x 2 matrix
    % "centroids"
    centroids = [centroidRow,centroidCol];
end