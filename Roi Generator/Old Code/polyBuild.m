function polyBuild(raw,pixlength,xres,filnm,innum1,innum2,innum3,innum4,TarDir,printQual)
%Polygon Builder
%For converting binary images into simple polylines



%The first part of the code will follow a similar algorithm to the
%vesselmosaic function we originally used to make quadrilaterals. Running
%this part of the code will attempt to find "points of interest" by first
%segmenting the image into 2 interleaved images, and then iteratively
%building quadrilaterals. This process will repeat for different values of
%quadrilateral heights, and points that appear more than once will be kept.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Begin editted code from Vesselmosaic
currdirec=cd;
dirr=cd;
a=raw;

temp=size(a);
ysize=temp(1);
xsize=temp(2);
b=a;

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

save('polyBuildOut.mat');

%End editted code from vesselmosaic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%