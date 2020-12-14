% Input variable images is a 3D image stack whose dimensions correspond to
% (rows,columns,z-slice)
function [filtMasks,ppOptions] = preprocess(images,metadata,auto) %,centroids
disp('Processing Images')
if auto == 0
ppOptions = ppSelector();
else
ppOptions =  [{'none'}    {[7]}    {[11]}]  ;
end
if ismember(1,ppOptions{1})==1
    images = uint16(removeLarge2(images,(metadata.scalingX*1000000)/0.1625));
    rawFile = [metadata.filename,' Raw.tif']; %metadata.filepath,
    if exist(rawFile,'file')
        delete(rawFile)
    end
    
    for i = 1:size(images,3)
        thisImg = uint16(images(:,:,i));
        imwrite(thisImg,rawFile,'WriteMode','append');
    end
else
    images = imresize(images,(metadata.scalingX*1000000)/0.1625);
    rawFile = [metadata.filename,' Raw.tif']; %metadata.filepath,
    if exist(rawFile,'file')
        delete(rawFile)
    end
    
    for i = 1:size(images,3)
        thisImg = uint16(images(:,:,i));
        imwrite(thisImg,rawFile,'WriteMode','append');
    end
end

filtMasks = (bpass3dMB(images, [1 1 1], [ppOptions{2} ppOptions{2} ppOptions{3}],[0 0]));
disp('Done Processing Images')

end