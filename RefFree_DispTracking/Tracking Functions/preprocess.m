% Input variable images is a 3D image stack whose dimensions correspond to
% (rows,columns,z-slice)
function [filtMasks,ppOptions] = preprocess(images,metadata) %,centroids
ppOptions = ppSelector();
if ismember(1,ppOptions{1})==1
    images = uint16(removeLarge2(images,(metadata.scalingX*1000000)/0.1625));
    rawFile = [metadata.filepath,metadata.filename,'Raw2.tif'];
    if exist(rawFile,'file')
        delete(rawFile)
    end
    
    for i = 1:size(images,3)
        thisImg = uint16(images(:,:,i));
        imwrite(thisImg,rawFile,'WriteMode','append');
    end
else
    images = imresize(images,(metadata.scalingX*1000000)/0.1625);
    rawFile = [metadata.filepath,metadata.filename,'Raw2.tif'];
    if exist(rawFile,'file')
        delete(rawFile)
    end
    
    for i = 1:size(images,3)
        thisImg = uint16(images(:,:,i));
        imwrite(thisImg,rawFile,'WriteMode','append');
    end
end
filtMasks = (bpass3dMB(images, [1 1 1], [ppOptions{2} ppOptions{2} ppOptions{3}],[0 0]));

end