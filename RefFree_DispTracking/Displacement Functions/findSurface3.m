function [Surface2,SurfaceAll] = findSurface3(image,raw,plane)
%This function attempts to build a surface by dividing the orginal image
%stack into simpler grids and determining when intensity falls below a
%particular threshold in Z. Only non-deformed regions are used to find the 
%surface, and the original image stack must have 'empty' frames at the 
%surface for this to work.

        %%
        clear maxStack
        filtS = 9;
        SE = strel('disk',filtS,6);        
        for i = 1:size(image.MaskStack,3)            
            sumInt(i) = sum(sum(image.MaskStack(:,:,i))); %find the total intensity per slice (to find minimum later)            
            maxStack(:,:,i) = ordfilt2(image.MaskStack(:,:,i),nnz(SE.Neighborhood),SE.Neighborhood); %apply a 2D max filter to the stack
            tStack(:,:,i) = maxStack(:,:,i)>10; %create a 3d thresholded version of the filtered stack
        end

        %%
        emptyslice = find(sumInt == min(sumInt)); %The dimmest slice is assumed to only contain noise
        shellStack = bwperim(tStack,26); %create a thin border around the 3D mask
        shellStack(:,:,1) = 0; %removed border at the bottom
        shellStack(:,:,emptyslice:end) = 0; %remove any pixels above the 'empty' slice
                
        for i = 1:size(shellStack,3)
            shellStack(:,:,i) = shellStack(:,:,i) .* (image.Area>0); %remove pixels within the cell area
        end
        
        SE2 = strel('disk',15,6); 
        borders(:,:) = imerode(shellStack(:,:,2)==0,SE2); %use this to delete border pixels in all images
        
        for i = 1:size(shellStack,3)
            shellStack(:,:,i) = shellStack(:,:,i) .* borders;
        end      
              
        idx = find(shellStack>0); %grab indices of all remaining pixels
        [xi,yi,zi] = ind2sub(size(shellStack),idx); %convert to coordinates
        
        y = (max(xi)+1-xi)*raw.dataKey(9,1); %convert pixel units to microns
        x = yi*raw.dataKey(9,1);
        z = zi*raw.dataKey(10,1);
  
        %fit coordinates to a surface
        Surface2 = fit([x,y],z,'poly11');
        SurfaceAll = fit([x,y],z,'poly11');
    
   

end