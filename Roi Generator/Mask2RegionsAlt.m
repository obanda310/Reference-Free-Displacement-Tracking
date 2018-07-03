function [polygons2,polygons] = Mask2RegionsAlt(raw,pixlength,outputName,TarDir,finalRes,options)
%Converts a binary image into a .Regions file and a .ovl file to be read by
%the Zen application.


%Set variables based on input arguments
temp=size(raw);
ysize=temp(1);
xsize=temp(2);
% save('raw.mat','raw')

polygons = bwperimtrace_frontend(raw);

%% Dilate polygon vertices to the ZEN coordinate space
%Convert pixel coordinates into real coordinates in units of micrometers by
%scaling with pixlength
polygons2 = polygons;
for i = 1:size(polygons2,1)
    polygons2{i,1}(:,1) = (finalRes(1,1)/xsize)*pixlength.*polygons2{i,1}(:,1)*(10^-6);
    polygons2{i,1}(:,2) = (finalRes(1,2)/ysize)*pixlength.*polygons2{i,1}(:,2)*(10^-6);
end


if (options(7,1) == 0) || (options(7,1) == 1 && mod(options(8,1),10)==1)
%% View all polylines as Raw Image Overlay
hFig = figure;
set(hFig, 'Position', [50 50 1200 600]);
subplot(1,2,1);
rawImage = max(raw,[],3);
imshow(rawImage);
hold on
for i = 1:size(polygons,1)
    try
        plot([polygons{i,1}(:,1);polygons{i,1}(1,1)],[polygons{i,1}(:,2);polygons{i,1}(1,2)])
        scatter([polygons{i,1}(:,1);polygons{i,1}(1,1)],[polygons{i,1}(:,2);polygons{i,1}(1,2)])
    end
end

%% View all polylines in Zen Frame of Reference
hFigS = subplot(1,2,2);
set(hFigS,'Position',[.55 .1 .4 .8]);
axis([-finalRes(1,1)/2*pixlength, finalRes(1,1)/2*pixlength, -finalRes(1,1)/2*pixlength, finalRes(1,1)/2*pixlength]);
xlabel('Microns')
ylabel('Microns')
hold on
for i = 1:size(polygons2,1)
    try
        plot([polygons2{i,1}(:,1);polygons2{i,1}(1,1)]/(10^-6),[polygons2{i,1}(:,2);polygons2{i,1}(1,2)]/(-1*10^-6))
        scatter([polygons2{i,1}(:,1);polygons2{i,1}(1,1)]/(10^-6),[polygons2{i,1}(:,2);polygons2{i,1}(1,2)]/(-1*10^-6))
    end
end
for i = 1:size(polygons2,1)
    xmin(i,1) = min(polygons2{i,1}(:,1));
    xmax(i,1) = max(polygons2{i,1}(:,1));
    ymin(i,1) = min(polygons2{i,1}(:,2));
    ymax(i,1) = max(polygons2{i,1}(:,2));
end
xMin = min(xmin);
xMax = max(xmax);
yMin = min(ymin);
yMax = max(ymax);
xwidth = (xMax-xMin) *10^6;
ywidth = (yMax-yMin) *10^6;

text(-finalRes(1,1)/2*pixlength+5,-finalRes(1,2)/2*pixlength+25,['Width: ',num2str(xwidth)])
text(-finalRes(1,1)/2*pixlength+5,-finalRes(1,2)/2*pixlength+10,['Height: ',num2str(ywidth)])
end


%% Create the .Regions file for ZEN
%Define path and file name for .rls file
filnm2=strcat(TarDir,'\',outputName,'.Regions');
filnm3=strcat(TarDir,'\',outputName,'.ovl');

%Pass the coordinates for the vertices of each ROI to the subroutine
%Poly2Regions or Poly2Overlay which will write the information about the
%list of ROIs to a file written in a format that is readable by Zeiss AIM
%or Zen software
Poly2Regions(filnm2,polygons2);
Poly2OVL(filnm3,polygons2);

%%
disp(['Regions file #' num2str(options(8,1)) ' has ',num2str(size(polygons,1)),' Regions and took ' num2str(toc) ' seconds!'])
%End of program


end