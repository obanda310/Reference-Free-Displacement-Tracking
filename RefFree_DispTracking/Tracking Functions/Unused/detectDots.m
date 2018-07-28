% Feature detection
function [boundaries,centroids] = detectDots(img)
    % find edges using Canny filter (log filter almost works too)
%     [bw,~] = edge(img,'log');
    [bw,~] = edge(img,'Canny',[0.1375,0.3]);
    % fill in gaps in detected edges by dilating features slightly
    se90 = strel('line',1,90);
    se0 = strel('line',1,0);
    % pack the image to make dilation operation more efficient - not 
    % necessary, but could speed up operations on large/many images
    % see: http://www.mathworks.com/help/images/ref/bwpack.html
    bwPacked = bwpack(bw);
    bwPackedDilate = imdilate(bwPacked,[se90,se0],'ispacked');
    bwDilate = bwunpack(bwPackedDilate);
    % fill in edge-detected areas so they are solid objects
    bwFill = imfill(bwDilate,'holes');
    % trace boundary of each object
    [boundaries,~] = bwboundaries(bwFill);
    % get centroid location of each object (in pixel units) 
    c = regionprops(bwFill,'Centroid');
    % concatenate all centroids into single matrix
    centroids = cat(1,c.Centroid);
end