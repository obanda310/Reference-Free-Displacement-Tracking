function output = vesselmosaic3v2(raw,pixlength,xres,filnm,innum1,innum2,innum3,innum4,TarDir,printQual)
% Optimized by Omar A. Banda 03/05/2018
% Major changes:
% -removed input option dialogues
% -replaced segment.m function with matlab's BWLABEL
% -removed generation of diagnostic tif image
% -

%Original Version Information:
% vesselmosaic3.m
%
% Copyright 2009-2011,
% Baylor College of Medicine and Rice University,
% All Rights Reserved
%
% ***Main program; run at the command prompt
%
% Reconstructs a binary image using a mosaic of quadrilateral ROIs to
% guide LSM photolithography.
%
% Two files are created in the working directory. The first file
% JCCxxxx.rls is an overlay file that defines a list of regions of
% interest (ROIs) in a format that is readable by Zeiss AIM or Zen
% software. By default, the overlay file is created with a .rls
% extension, which is the extension used by the AIM software. For use
% with the Zen software, simply change the file extension to .ovl
%
% The second file ROIxxxx.tif is a rendering of the ROIs in tif format.
% This is for making qualitative comparisons between the generated ROIs
% and the original image, but is not intended for quantitative
% comparisons.
%
% .rls files must be placed into the folder C:\AIM\ROIS
% .ovl files must be loaded using the Zen GUI
%
% Input arguments:
% raw = Matrix of integers, all either 0 or 1, describing the
% binary image to be used for image guided patterning
% pixlength = Physical length of a single pixel, in micrometers
% xres = Desired maximum width of an ROI, in pixels. Must be an
% integer > 1
% filnm = Name and path of a template .rls file. This file must be
% a .rls file defining a list of Polyline ROIs
% innum1 = First digit in file identifier number
% innum2 = Second digit in file identifier number
% innum3 = Third digit in file identifier number
% innum4 = Fourth digit in file identifier number
%
% Notes:
% - file identifier number is an arbitrary tag for distinguising between
% different overlay files
%
%Original Version Information (end)

%Set variables based on input arguments
strinnum1=int2str(innum1);
strinnum2=int2str(innum2);
strinnum3=int2str(innum3);
strinnum4=int2str(innum4);

currdirec=cd;
dirr=cd;

a=raw;

%Display original image for reference
figure('Position',[1 1 659 659]);
imagesc(a);colormap(gray);
set(gca,'DataAspectRatio',[1 1 1]);

temp=size(a);
ysize=temp(1);
xsize=temp(2);
b=a;

%% OAB: We never use these options, removing for now
% button='';
% button = questdlg('Shrink or Expand?','Scaling','Yes','No','No');
% [dummy1 dummy2]=size(button);
%
% if dummy2==3
%     indata=inputdlg({'Area Scaling (% of original)'},'Scaling parameters',1,{'100'});
%     relsize=str2num(indata{1});
%     if relsize~=100
%         relsize=relsize/100;
%
%         map=segment(b);
%         map2=map(:,:,1).*b;
%         map3=map2.*0;
%         for i=1:(max(max(map2)))
%             [dummyrow dummycol]=size(map2);
%             curr=map2==i;
%             origpix=sum(sum(curr));
%             currpix=origpix;
%             if origpix==0
%                 continue
%             end
%             if relsize<1
%                 while (currpix>(relsize*origpix)) && (currpix~=0)
%                     curr=bwmorph(curr,'erode',1);
%                     currpix=sum(sum(curr));
%                 end
%             elseif relsize>1
%                 while (currpix<(relsize*origpix)) && (currpix~=(dummyrow*dummycol))
%                     curr=bwmorph(curr,'dilate',1);
%                     currpix=sum(sum(curr));
%                 end
%             end
%             map3=map3+curr;
%         end
%         map3=map3>0;
%         a=map3;
%         b=a;
%     end
% end
%%

%button='';
%button = questdlg('Select ROI Orientation','Orientation','Vertical','Horizontal','Vertical');
%[dummy1 dummy2]=size(button);
dummy2 = 10; %we never use anything but Horizontal
if dummy2==10
    
    %HORIZONTAL%
    
    %Divide input image into two interleaved images before segmentation. This
    %essentially divides large regions into smaller, simpler regions
    for i=1:(2*xres):ysize
        b(i:(i+xres-1),:)=0;
    end
    b=b(1:ysize,1:xsize);
    if (ysize-i)>xres
        b((ysize-xres+1):ysize,:)=0;
    end
    c=(1-b);
    d=c.*a;
    
    
    %Pass each of the two interleaved images to the subroutine segment.m to
    %identify and count each interconnected region
    
    %e=segment(b);
    e = bwlabel(b);
    %f=segment(d);
    f = bwlabel(d);
    g=e(:,:,1).*b;
    h=f(:,:,1).*d;
    
    
    %Initialize matrix of coordinates
    splinex=[0 0 0 0];
    spliney=[0 0 0 0];
    splinenum=1;
    
    %Approximate each segmented region into a vertically simple quadrilateral.
    %Begin with the first of the two interleaved images.
    
    %Display cumulative progress as a percentage
    
    for i=1:(max(max(g)))
        DefiningROIs=100*i/(max(max(g)));
        if sum(sum(g==i))==0
            continue
        end
        curr=g==i;
        [row col]=find(curr==1);
        
        %Define upper and lower limits of quadrilateral
        uy=min(row);
        ly=max(row);
        
        %Define left and right limits of quadrilateral at each y-limit
        [row1 col1]=find(curr(uy,:)==1,1,'first');
        [row2 col2]=find(curr(uy,:)==1,1,'last');
        ulx=col1;
        urx=col2;
        if ulx==urx
            urx=ulx+1;
        end
        
        [row1 col1]=find(curr(ly,:)==1,1,'first');
        [row2 col2]=find(curr(ly,:)==1,1,'last');
        llx=col1;
        lrx=col2;
        if llx==lrx
            lrx=llx+1;
        end
        
        %Store the coordinates of these vertices into splinex and spliney
        splinex(splinenum,:)=[ulx urx lrx llx];
        spliney(splinenum,:)=[uy uy ly ly];
        splinenum=splinenum+1;
        
    end
    
    
    %Now, repeat the above process for the second of the two interleaved images
    for i=1:(max(max(h)))
        DefiningROIs2=100*i/(max(max(h)));
        if sum(sum(h==i))==0
            continue
        end
        curr=h==i;
        [row,col]=find(curr==1);
        uy=min(row);
        ly=max(row);
        
        [~,col1]=find(curr(uy,:)==1,1,'first');
        [~,col2]=find(curr(uy,:)==1,1,'last');
        
        ulx=col1;
        urx=col2;
        if ulx==urx
            urx=ulx+1;
        end
        
        [~,col1]=find(curr(ly,:)==1,1,'first');
        [~,col2]=find(curr(ly,:)==1,1,'last');
        llx=col1;
        lrx=col2;
        if llx==lrx
            lrx=llx+1;
        end
        splinex(splinenum,:)=[ulx urx lrx llx];
        spliney(splinenum,:)=[uy uy ly ly];
        splinenum=splinenum+1;
        
    end
    
else
    %%
    % %VERTICAL%
    %
    % %Divide input image into two interleaved images before segmentation. This
    % %essentially divides large regions into smaller, simpler regions
    % for i=1:(2*xres):xsize
    % b(:,i:(i+xres-1))=0;
    % end
    % b=b(1:ysize,1:xsize);
    % if (xsize-i)>xres
    % b(:,(xsize-xres+1):xsize)=0;
    % end
    % c=(1-b);
    % d=c.*a;
    %
    % %Pass each of the two interleaved images to the subroutine segment.m to
    % %identify and count each interconnected region
    % e=segment(b);
    % f=segment(d);
    % g=e(:,:,1).*b;
    % h=f(:,:,1).*d;
    %
    % %Initialize matrix of coordinates
    % splinex=[0 0 0 0];
    % spliney=[0 0 0 0];
    % splinenum=1;
    %
    % %Approximate each segmented region into a vertically simple quadrilateral.
    % %Begin with the first of the two interleaved images.
    %
    % %Display cumulative progress as a percentage
    % for i=1:(max(max(g)))
    %     DefiningROIs=100*i/(max(max(g)))
    % if sum(sum(g==i))==0
    % continue
    % end
    % curr=g==i;
    % [row col]=find(curr==1);
    %
    % %Define left and right limits of quadrilateral
    % lx=min(col);
    % rx=max(col);
    %
    % %Define upper and lower limits of quadrilateral at each x-limit
    % [row col]=find(curr(:,lx)==1);
    % uly=min(row);
    % lly=max(row);
    % if uly==lly
    % lly=lly+1;
    % end
    % [row col]=find(curr(:,rx)==1);
    % ury=min(row);
    % lry=max(row);
    % if ury==lry
    % lry=lry+1;
    % end
    %
    % %Store the coordinates of these vertices into splinex and spliney
    % splinex(splinenum,:)=[lx rx rx lx];
    % spliney(splinenum,:)=[uly ury lry lly];
    % splinenum=splinenum+1;
    %
    % %Go to the next region and repeat
    % end
    %
    % %Now, repeat the above process for the second of the two interleaved images
    % for i=1:(max(max(h)))
    %     DefiningROIs2=100*i/(max(max(h)))
    % if sum(sum(h==i))==0
    % continue
    % end
    % curr=h==i;
    % [row col]=find(curr==1);
    % lx=min(col);
    % rx=max(col);
    % [row col]=find(curr(:,lx)==1);
    % uly=min(row);
    % lly=max(row);
    % if uly==lly
    % lly=lly+1;
    % end
    % [row col]=find(curr(:,rx)==1);
    % ury=min(row);
    % lry=max(row);
    % if ury==lry
    % lry=lry+1;
    % end
    % splinex(splinenum,:)=[lx rx rx lx];
    % spliney(splinenum,:)=[uly ury lry lly];
    % splinenum=splinenum+1;
    % end
    
end
%%
%Here, the program pauses - to continue, type "return" at the command
%prompt
%keyboard

% OAB removed to save time
%MATLAB is used to create a rendering of the ROIs. This is for making
%qualitative comparisons between the generated ROIs and the original image,
%but is not intended for quantitative comparisons.
if printQual == 1
j=a.*0;
temp=size(splinex);
count = 0;
for i=1:temp(1)
    count = count + 1;
    PlottingROIs=100*i/temp(1);
    BW=roipoly(a,(splinex(i,:)+[-1 0 0 -1]),(spliney(i,:)+[-1 -1 0 0]));
    j=j+BW;
end
j=j>0;
figure('Position',[1 1 659 659]);
imagesc(j);colormap(gray);
set(gca,'DataAspectRatio',[1 1 1])

%Rendering is saved in the working directory as 'ROIxxxx.tif'
temp=strcat(TarDir,'\','ROI',strinnum1,strinnum2,strinnum3,strinnum4,'.tif');
imwrite(j,temp,'tif');
end
%%

%Convert pixel coordinates into real coordinates in units of micrometers by
%scaling with pixlength
splinextrans=pixlength.*splinex;
splineytrans=pixlength.*spliney;

%Convert into the coordinate system used by the microscope. Coordinate
%system is in units of meters, with the origin at the center of the field
%of view.
splinextrans=((splinextrans-(xsize*pixlength/2))*(10^-6));
splineytrans=((splineytrans-(ysize*pixlength/2))*(10^-6));

%Define path and file name for .rls file
filnm2=strcat(TarDir,'\','AS',strinnum1,strinnum2,strinnum3,strinnum4,'.Regions');
tempsizee=size(splinextrans);
drawnumber=tempsizee(1);

%Pass the coordinates for the vertices of each ROI to the subroutine
%roifiles2.m which will write the information about the list of ROIs to a
%file written in a format that is readable by Zeiss AIM or Zen software
regionfiles2(filnm,filnm2,strinnum1,strinnum2,strinnum3,strinnum4,drawnumber,splinextrans,splineytrans);
cd(currdirec);

%Some of the important variables are included in the output for further use
%at the command prompt, if necessary.
output.original=a;
%output.roi=j;
output.splinex=splinex;
output.spliney=spliney;
output.splinextrans=splinextrans;
output.splineytrans=splineytrans;
%output.count=count;
disp(['Mask Has ',num2str(drawnumber),' Regions and took ' num2str(toc) ' seconds!']) 
%End of program
