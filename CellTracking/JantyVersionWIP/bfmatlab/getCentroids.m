function centroids = getCentroids(img,mask)
    % Use regionprops to find the weighted (i.e. by grayscale intensity in the original image)
    % centroids of all the objects in the new filtered masks
    % Use regionprops to find the weighted centroids of all the objects in
    % the image. The objects are determined using the mask that corresponds
    % to the image, and the weighted locations of the centroids are
    % determined using the grayscale pixel intensity values (obtained from
    % the original image) within each object
    c = regionprops(mask,img,'WeightedCentroid');
    % Format the centroid information so that each cell in "centroids"
    % is a 2-column matrix that has the x and y positions for each
    % centroid in the image that corresponds to that cell
    centroids = cat(1,c.WeightedCentroid);
end