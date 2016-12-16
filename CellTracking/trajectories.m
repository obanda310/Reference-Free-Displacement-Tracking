close all; clear;
%Analyzing Trajectories from FIJI input or from Custom Code

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Section 1.) Inputs and Sorting the Raw Data From Excel File
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1.1 Loading Data and Background Images for Overlays%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[num,dataKey] = InputSelector();

[nameFluorFile,filePath] = uigetfile('*.tif','Select Fluorescent Image for Overlay');
imageFluor = imread([filePath,nameFluorFile]);

[nameTransFile,filePath] = uigetfile('*.tif','Select Transmitted Image for Overlay');
imageTrans = imread([filePath,nameTransFile]);

[nameBlackFile,filePath] = uigetfile('*.tif','Select a Black Image of the Correct Dimensions');
imageBlack = imread([filePath,nameBlackFile]);

outputs = OutputSelector();

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1.2 Input Variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numIndices = 6; %for Section 1 and 2 (number of elements in book1)
numTraj = max(num(:,dataKey(4,1))); %Number of Trajectories
totalNumFrames = max(num(:,dataKey(3,1)))+dataKey(8,1); %Maximum number of frames observable for any one object
book1 = zeros(numIndices,totalNumFrames,numTraj); %Creates book1
book2 = zeros(numTraj,9);
%%

for i = 1:numTraj
    % Here we build a book of pages (3D array) with the data for a single
    % object/trajectory per page. Most of the code is to ensure that each
    % page is the same size matrix as the next.
    
    % stores data to 'tempObj' pertaining to all frames of one
    % object/trajectory
    tempObj = num(num(:,dataKey(4,1))==i,:);
    % number of frames that the current object appears in.
    numFrames = size(tempObj,1);
    if numFrames > 0
        % the first frame that an object appears in (tracking software starts
        % at 0
        startFrame = min(tempObj(:,dataKey(3,1)))+dataKey(8,1);
        % the last frame that an object appears in
        endFrame = (startFrame+numFrames-1);
        book2(i,3) = startFrame;
        book2(i,4) = endFrame;
        % this fills in the upper portion of obj matrix with zeros if the first
        % frame is not 0
        
        book1(1,startFrame:endFrame,i) = tempObj(:,dataKey(1,1)).*dataKey(7,1);
        book1(2,startFrame:endFrame,i) = tempObj(:,dataKey(2,1)).*dataKey(7,1);
        book1(6,startFrame:endFrame,i) = tempObj(:,dataKey(5,1));
    else
        startFrame = 1;
        endFrame = 1;
        book2(i,3) = startFrame;
        book2(i,4) = endFrame;
    end
    if i == 1 || i == 5000 || i == 10000 || i == 15000 || i==20000 || i==25000
        disp(['Progress: ' num2str(i) ' of ' num2str(numTraj)])
    end
end
disp('done Section 1')
%Here we use the temporary data stored in 'obj' to build 2 matrices to
%describe the x and y positions in frame 1 for each object, as well as the
%displacement from those matrices occurring in successive frames.

%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Section 2.) Building book1 from the raw data
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


%Notes: The majority of indices are added in later sections, but listed
%here for reference.

%Index List for Book1 - Frame Dependent Values
%1 = Raw X Centroid Location of Traj in Specified Frame
%2 = Raw Y Centroid Location of Traj in Specified Frame
%3 = Raw dX Centroid Location of Traj in Specified Frame from First Frame
%4 = Raw dY Centroid Location of Traj in Specified Frame from First Frame
%5 = Magnitude of displacement
%6 = Intensity in current frame (column 14 from TrackMate output)

%Index List for Book2 - Frame Independent Values
%1 = Raw X Centroid Location of Traj in First Frame
%2 = Raw Y Centroid Location of Traj in First Frame
%3 = First Frame that a Traj Appears
%4 = Last Frame that a Traj Appears
%5 = Value of book1 index 5 (see above) in Last Frame that Traj Appears
%6 = Value of book1 index 6 (see above) in Last Frame that Traj Appears
%7 = Value of book1 index 1 (see above) in Last Frame that Traj Appears
%8 = Value of book1 index 2 (see above) in Last Frame that Traj Appears
%9 = Maximum Magnitude of displacement

for i = 1:numTraj
    book2(i,1) = book1(1,book2(i,3),i);
    book2(i,2) = book1(2,book2(i,3),i);
    book1(3,book2(i,3):book2(i,4),i) = book1(1,book2(i,3):book2(i,4),i) - book2(i,1);
    book1(4,book2(i,3):book2(i,4),i) = book1(2,book2(i,3):book2(i,4),i) - book2(i,2);
    book1(5,:,i) = (book1(3,:,i).^2 + book1(4,:,i).^2).^0.5;
    book2(i,10) = i;
end

%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Section 3.) Identifying Neighborhood Trajectories for Each Trajectory
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
clear tempDistances
book3 = zeros(numTraj,50);
cm3 = zeros(numTraj,50,totalNumFrames);
maxDistance = (max(max(book1(5,:,:))));
for i = 1:numTraj
    tempDistances(1:numTraj,1) = ((book2(:,1)-book2(i,1)).^2.+(book2(:,2)-book2(i,2)).^2).^0.5;
    tempDistances(1:numTraj,2) = linspace(1,numTraj,numTraj);
    tempDistances = sortrows(tempDistances);
    book3(i,1:50) = tempDistances(1:50,2);
    for j = 1:50
        if book2(book3(i,j),9) < (maxDistance/4)
            cm3(i,j,1:totalNumFrames) = book1(6,1:totalNumFrames,book3(i,j));
        else
            cm3(i,j,1:totalNumFrames) = NaN;
        end
    end
end

%%

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Section 4.)Identifying Deviations From Zero-State in Later Frames
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4.2 Storing the Maximum Displacement Values%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Notes: Some trajectories do not persist all the way to the top frame, and
%so they would not be visible in the outputs in section 5 unless their
%maximum value is used for those overlays.

for i = 1:numTraj
    
    book2(i,5) = book1(3,book2(i,4),i);
    book2(i,6) = book1(4,book2(i,4),i);
    book2(i,7) = book1(1,book2(i,4),i);
    book2(i,8) = book1(2,book2(i,4),i);
    book2(i,9) = (book2(i,5).^2 + book2(i,6).^2).^0.5;
end

disp('done Section 4.2')
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4.3 Creating a Color Map for Quiver Plots%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[cm1,cm2,cmD,cmDS,colorMap,colorScheme] = createColorMap(book1,book2);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4.4a Using Intensity Values to Extract Z-Information%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Notes: Calculates average distance between PSFs in a single pillar.
%Outputs stored as variable "yDataDiffAverage"
if ismember(10,outputs) == 1
    close all
    
    intProfiles = zeros(totalNumFrames,numTraj);
    
    sortedArray = zeros(numTraj,3);
    for i = 1:numTraj
        sortedArray(i,1) = i;
        sortedArray(i,2) = book2(i,1);
        sortedArray(i,3) = book2(i,2);
    end
    
    [sortedArray,sortedArrayIndices] = sortrows(sortedArray,[3, 2]);
    
    for i = 1:(numTraj)
        intProfiles(book2(sortedArray(i,1),3):book2(sortedArray(i,1),4),i) = book1(6,book2(sortedArray(i,1),3):book2(sortedArray(i,1),4),sortedArray(i,1));
    end
    %image(intProfiles,'cdatamapping','scaled');
    %colormap(gray);
    %imwrite(intProfilesImage,'intProfilesImage.tiff');
    close;
    intProfilesPS = phasesym(intProfiles,'minWaveLength',3,'mult',1.5);
    intProfilesMax = (localmax(intProfilesPS'))' > 0;
    intProfilesClosed1 = imclose(intProfilesMax, ones(1,20));
    intProfilesCleaned = bwareaopen(intProfilesClosed1,200);
    intProfilesClosed2 = imclose(intProfilesCleaned, ones(1,1000));
    
    intProfilesLabeled = bwlabel(intProfilesClosed2);
    
    intProfilesFitCoefs = zeros(5,5);
    yData = zeros(numTraj,5);
    yDataDiffAverage = zeros(1,5);
    intProfilesFits = figure;
    image(intProfiles,'cdatamapping','scaled'); colormap(gray);
    hold on
    for i = 1:max(max(intProfilesLabeled))
        [intProfilesCurrentLabely,intProfilesCurrentLabelx] = find((intProfilesLabeled == i).*intProfilesCleaned);
        intProfilesFitCoefs(1:2,i) = polyfit(intProfilesCurrentLabelx, intProfilesCurrentLabely, 1);
        xData = linspace(1,numTraj,numTraj);
        yData(1:numTraj,i) = intProfilesFitCoefs(1,i)*(xData) + intProfilesFitCoefs(2,i);
        plot(xData,yData(1:numTraj,i))
        hold on
        plot(intProfilesCurrentLabelx, intProfilesCurrentLabely,'*','MarkerSize',5);
        hold on
    end
    savefile = [filePath '\Intensity Profile Linear Fits.tif'];
    export_fig(intProfilesFits,savefile);
    hold off
    
    intProfilesRaw = figure;
    image(intProfiles,'cdatamapping','scaled'); colormap(gray);
    savefile = [filePath '\Intensity Profiles 1 Raw.tif'];
    export_fig(intProfilesRaw,savefile);
    
    intProfilesPSfig = figure;
    image(intProfilesPS,'cdatamapping','scaled'); colormap(gray);
    savefile = [filePath '\Intensity Profiles 2 Phase Sym.tif'];
    export_fig(intProfilesPSfig,savefile);
    
    intProfilesMaxfig = figure;
    image(intProfilesMax,'cdatamapping','scaled'); colormap(gray);
    savefile = [filePath '\Intensity Profiles 3 Max.tif'];
    export_fig(intProfilesMaxfig,savefile);
    
    intProfilesClosedfig1 = figure;
    image(intProfilesClosed1,'cdatamapping','scaled'); colormap(gray);
    savefile = [filePath '\Intensity Profiles 4 Closed.tif'];
    export_fig(intProfilesClosedfig1,savefile);
    
    intProfilesCleanedfig = figure;
    image(intProfilesCleaned,'cdatamapping','scaled'); colormap(gray);
    savefile = [filePath '\Intensity Profiles 5 Cleaned.tif'];
    export_fig(intProfilesCleanedfig,savefile);
    
    intProfilesClosedfig2 = figure;
    image(intProfilesClosed2,'cdatamapping','scaled'); colormap(gray);
    savefile = [filePath '\Intensity Profiles 6 Closed.tif'];
    export_fig(intProfilesClosedfig2,savefile);
    
    intProfilesLabeledfig = figure;
    image(intProfilesLabeled,'cdatamapping','scaled'); colormap(gray);
    savefile = [filePath '\Intensity Profiles 7 Labeled.tif'];
    export_fig(intProfilesLabeledfig,savefile);
    
    for i = 2:max(max(intProfilesLabeled))
        yDataDiffAverage(1,i) = mean(yData(:,i-1)-yData(:,i));
    end
    
    
    
    
    
    % close;
    % intProfilesPSM = phasesymmono(intProfiles,'minWaveLength',2,'mult',1.5);
    % intProfilesPS = phasesym(intProfiles,'minWaveLength',2,'mult',1.5);
    % intProfilesPCM = phasecongmono(intProfilesPSM);
    % [intProfilesPC3,trash, orient] = phasecong3(intProfilesPS);
    % [intProfilesNMS, intProfilesNMSPos] = nonmaxsup(intProfilesPS,intProfilesOrient,1.3);
    % intProfilesNMS2 = medfilt2(intProfilesNMS, [3 3]) > 0;
    % intProfilesNMS3 = imclose(intProfilesNMS2, ones(1,200));
    % intProfilesNMSM = nonmaxsuppts(intProfilesPCM);
    %
    % image(intProfilesNMS3,'cdatamapping','scaled'); colormap(gray);
    
    %imwrite(intProfilesImage,'intProfilesImage.tiff');
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4.4b Using Intensity Values to Extract Z-Information%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(9,outputs) == 1
    cm2(cm2 == 0) = NaN;
    % plottedProfiles = figure
    % for i = 1:numTraj
    % plot(linspace(1,54,54),cMDBook2(30,:,i,18),'.','MarkerSize',5)
    % hold on
    % end
    
    trajOfInterest = 2500;
    %     winSize = 50;
    %     if trajOfInterest > winSize && trajOfInterest < (numTraj-winSize)
    %         lowerLim = trajOfInterest - winSize;
    %         upperLim = trajOfInterest + winSize;
    %     elseif trajOfInterest < winSize
    %         lowerLim = trajOfInterest-(trajOfInterest-1);
    %         upperLim = trajOfInterest*2 + winSize;
    %     elseif trajOfInterest > (numTraj-winSize)
    %         lowerLim = trajOfInterest - winSize - (numTraj-trajOfInterest);
    %         upperLim = trajOfInterest + (trajOfInterest-1);
    %     end
    
    plottedProfiles = figure
    subplot(2,1,1)
    plot(linspace(1,totalNumFrames,totalNumFrames),book1(6,:,trajOfInterest),'.','MarkerSize',5,'Color',[0 0 0])
    hold on
    findpeaks((book1(6,:,trajOfInterest)))
    
    currentAverages = zeros(totalNumFrames,1);
    
    for i = 1:totalNumFrames
        errorbar(i,mean(cm3(trajOfInterest,1:50,i),'omitnan'),std(cm3(trajOfInterest,1:50,i),'omitnan'),'.','MarkerSize',5,'Color','r');
        currentAverages(i,1) = mean(cm3(trajOfInterest,1:50,i),'omitnan');
    end
    currentAverages(isnan(currentAverages))=0;
    findpeaks(currentAverages(:,1))
    
    hold on
    xVals = linspace(1,totalNumFrames,totalNumFrames);
    fitWeights = zeros(totalNumFrames,1);
    fitWeights(:,1) = xVals(1,:)+1000;
    s = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[4000,25000,0,0,0,0,0],...
        'Upper',[6000,35000,15,.1,1,3.14,15],...
        'StartPoint',[5500,29000,11.5,.015,0.1731,2.2,6.712]);
    f = fittype('pillarFit2(x,S,iS,iP,P,Pg,Mo,Mf)','options',s);
    [fitCoef,fitQual] = fit(xVals',currentAverages,f,'weight',fitWeights);
    fitCoef2 = coeffvalues(fitCoef);
    yVals = pillarFit(xVals,fitCoef2(1,5),fitCoef2(1,7),fitCoef2(1,6),fitCoef2(1,3),fitCoef2(1,4),fitCoef2(1,2),fitCoef2(1,1));
    plot(xVals,yVals,'MarkerSize',20,'Color','bl')
    findpeaks(yVals(1,:))
    
    hold off
    
    subplot(2,1,2)
    imshow(imageBlack)
    hold on
    for i = 1:size(book3,2)
        scatter(book2(book3(trajOfInterest,i),1),book2(book3(trajOfInterest,i),2),'.','g')
    end
    scatter(book2(trajOfInterest,1),book2(trajOfInterest,2),'o','b')
    
    
    % folderName = strcat( 'Intensity Average Profiles with ',num2str(cmD),' Divisions');
    % mkdir(blackPath,folderName)
    % for j = 1:cmD
    % plottedProfilesErrorBars = figure
    % for i = 1:totalNumFrames
    % errorbar(i,mean(cMDBook2(30,i,:,j),'omitnan'),std(cMDBook2(30,i,:,j),'omitnan'),'.','MarkerSize',5)
    % hold on
    % end
    %
    % axis([0 60 0 40000])
    % savefile = [blackPath '\' folderName '\Displacement Magnitude Division ' num2str(j) '.tif'];
    % export_fig(plottedProfilesErrorBars,savefile);
    %
    % close
    % end
end
%%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%5.)IMAGE OUTPUTS
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%5.2 Drawing Zero-State Displacement Fields on Black Background%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ismember(2,outputs) == 1
    folderName = strcat('Black Image_',colorScheme,' Quiver Overlays');
    mkdir(filePath,folderName)
    for f = 1:totalNumFrames        % number of z-slices
        blackOverlay = figure('Position',[0 0 1000 1000]);
        imshow(imageBlack,[])
        hold on
        
        for i = 1:cmD
            quiver(cm1(1,f,:,i),cm1(2,f,:,i),cm1(3,f,:,i),cm1(4,f,:,i),...
                0,'color',[colorMap(i,1:3)]);
            hold on
        end
        hold off
        savefile = [filePath '\' folderName '\Black Background Overlay' colorScheme ' ' num2str(f) '.tif'];
        export_fig(blackOverlay,savefile,'-native');
        close
    end
end



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%5.3 Plotting Centroids on a Black Background%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Notes: This section plots centroids (as in 5.3) on a black background for
%processing into heat maps in imageJ to view 3D information in a 2D image.

if ismember(4,outputs) == 1
    folderName = strcat( 'Centroid_',colorScheme,' Black Image Overlays');
    mkdir(filePath, folderName)
    for f = 1:totalNumFrames
        centroidsOnly = figure('Position',[0 0 1000 1000]);
        imshow(imageBlack,[])
        hold on
        for i = 1:cmD
            
            xTemp = squeeze(squeeze(cm1(1,f,:,i)));
            xTemp = xTemp';
            yTemp = squeeze(squeeze(cm1(2,f,:,i)));
            yTemp = yTemp';
            plot(xTemp,yTemp,'.','MarkerSize',3,'Color',[colorMap(i,1:3)]);
            
            hold on
        end
        hold on
        hold off
        savefile = [filePath '\' folderName '\Centroids on Frame ' colorScheme ' ' num2str(f) '.tif'];
        export_fig(centroidsOnly,savefile,'-native');
        close
    end
end

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%5.4 DEBUGGING:Plotting Corrected X,Y Coordinate Fields%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: This sections is really helpful for debugging. It shows x0,y0 positions
% as well as current X,Y coordinates for each trajectory on each frame. It
% also plots the final displacement quiver field.

if ismember(3,outputs) == 1
    debugImage = figure('Position',[0 0 1000 1000]);
    imshow(imageTrans,[])
    hold on
    xTemp = squeeze(book1(1,1,:));
    yTemp = squeeze(book1(2,1,:));
    plot(xTemp,yTemp,'r.','MarkerSize',10);
    hold on
    if ismember(7,outputs) == 0
        for f = 1:totalNumFrames
            xTemp = squeeze(book1(1,f,:));
            yTemp = squeeze(book1(2,f,:));
            plot(xTemp,yTemp,'b.','MarkerSize',3);
            hold on
        end
    end
    hold on
    if ismember(8,outputs) == 0
        for i = 1:cmD
            quiver(cm2(:,1,i),cm2(:,2,i),cm2(:,5,i),cm2(:,6,i),0,'color',[colorMap(i,1:3)]);
            hold on
        end
        
    end
    hold off
    savefile = [filePath '\Debug Image ' colorScheme '.tif'];
    export_fig(debugImage,savefile,'-native');
end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%5.5 Plotting Transmitted Quiver Overlays v1%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Notes: 5.5 Saves files, 5.6 does not.

if ismember(5,outputs) == 1
    trajOverlay = figure('Position',[0 0 1000 1000]);
    imshow(imageFluor,[])
    hold on
    for i = 1:cmD
        
        quiver(cm2(:,1,i),cm2(:,2,i),cm2(:,5,i),cm2(:,6,i),0,'color',[colorMap(i,1:3)]);%(i^2)/((cmD/1.5)^2) also (i^1.3)/((cmD/2)^1.3)
        hold on
    end
    %quiver(book1(10,1,:),book1(11,1,:),book1(12,totalNumFrames,:),book1(13,totalNumFrames,:),0,'g');
    hold off
    savefile = [filePath '\Fluorescent Overlay ' colorScheme '.tif'];
    export_fig(trajOverlay,savefile,'-native');
end

%%
if ismember(5,outputs) == 1
    transmittedOverlay = figure('Position',[0 0 1000 1000]);
    imshow(imageTrans,[])
    hold on
    for i = 1:cmD
        quiver(cm2(:,1,i),cm2(:,2,i),cm2(:,5,i),cm2(:,6,i),0,'color',[colorMap(i,1:3)]);%(i^1.3)/((cmD/2)^1.3)
        hold on
    end
    %quiver(book1(10,1,:),book1(11,1,:),book1(12,totalNumFrames,:),book1(13,totalNumFrames,:),0,'g');
    hold off
    savefile = [filePath '\Transmitted Overlay ' colorScheme '.tif'];
    export_fig(transmittedOverlay,savefile,'-native');
end

%%
clear MeshBook MeshList
clear book4

book4 = book1;
book4(~book4) = NaN;

MeshBook = book4(1:2,1:25,1:121);
MeshBook(2,:,:) = MeshBook(2,:,:)*-1;
node = 0;

for i = 1:size(MeshBook,2)
    for j = 1:size(MeshBook,3)
        node = node + 1;
        if isnan(MeshBook(1,i,j)) == 1           
            MeshBook(1,i,j) = MeshBook(1,i-1,j);      
        end
        if isnan(MeshBook(2,i,j)) == 1                        
            MeshBook(2,i,j) = MeshBook(2,i-1,j);
        end
        MeshBook(3,i,j) = i;
        MeshBook(4,i,j) = node;
        
        
        MeshList(2,node) = MeshBook(1,i,j); 
        MeshList(3,node) = MeshBook(2,i,j);
        MeshList(4,node) = MeshBook(3,i,j);
        MeshList(1,node) = MeshBook(4,i,j);
        
        MeshList2(2,node) = MeshBook(1,1,j); 
        MeshList2(3,node) = MeshBook(2,1,j);
        MeshList2(4,node) = MeshBook(3,i,j);
        MeshList2(1,node) = MeshBook(4,i,j);
    end
end

%Building Elements
%Start with the first frame
clear EleList EleList2
    Element = 0;
for i = 1:size(MeshBook,3)
    clear current diff
    current(1:2,:) = MeshBook(1:2,1,:);
    diff(1,:) = current(1,:) - MeshBook(1,1,i);
    diff(2,:) = current(2,:) - MeshBook(2,1,i);
    diff(3,:) = sqrt(diff(1,:).^2.+diff(2,:).^2);
    diff(diff<-5) = NaN;
    diff(diff>22) = NaN;
    for j = 1:size(diff,2)
        if isnan(diff(1,j))==1 || isnan(diff(2,j))== 1 || isnan(diff(3,j)) ==1 || diff(3,j) == 0
            diff(1:3,j) = 0;
        end
    end
    
    clear current2
    current2 = find(diff(3,:));
    
    if size(current2,2) == 3
        current3 = diff(1:3,current2);
         mesh3 = current2(find(current3(3,:)==max(current3(3,:))));
         mesh4 = current2(find(current3(2,:)==min(current3(2,:))));
         mesh2 = current2(find(current3(1,:)==min(current3(1,:))));
        for k = 1:(size(MeshBook,2)-1)
        Element = Element + 1;
        
        EleList(i,1,k) = MeshBook(4,k,i); %5
        EleList(i,2,k) = MeshBook(4,k,mesh2); %6
        EleList(i,3,k) = MeshBook(4,k,mesh3); %7
        EleList(i,4,k) = MeshBook(4,k,mesh4); %8
        EleList(i,5,k) = MeshBook(4,k+1,i); %1
        EleList(i,6,k) = MeshBook(4,k+1,mesh2); %2
        EleList(i,7,k) = MeshBook(4,k+1,mesh3); %3
        EleList(i,8,k) = MeshBook(4,k+1,mesh4); %4

        EleList(i,9,k) = Element;
        EleList2(2:9,Element) = EleList(i,1:8,k);
        EleList2(1,Element) = Element;
        end
    end
        
    
end
%%

meshTxt = fopen('mesh.txt','wt');
nodesFormat = '<node id=" %d "> %f, %f, %f </node>\n';
fprintf(meshTxt,'<?xml version="1.0" encoding="ISO-8859-1"?>\n<febio_spec version="2.5">\n<Geometry>\n<Nodes name="Part1">\n');
fprintf(meshTxt,nodesFormat,MeshList(1:4,:));
fprintf(meshTxt,'</Nodes>\n<Elements type="hex8" mat="1" name="part1">\n');
elementsFormat = '<elem id=" %d "> %d, %d, %d, %d, %d, %d, %d, %d </elem>\n';
fprintf(meshTxt,elementsFormat,EleList2(1:9,:));
fprintf(meshTxt,'</Elements>\n</Geometry></febio_spec>');
fclose(meshTxt);
%%
%Undeformed Mesh

meshTxt = fopen('meshUD.txt','wt');
nodesFormat = '<node id=" %d "> %f, %f, %f </node>\n';
fprintf(meshTxt,'<?xml version="1.0" encoding="ISO-8859-1"?>\n<febio_spec version="2.5">\n<Geometry>\n<Nodes name="Part1">\n');
fprintf(meshTxt,nodesFormat,MeshList2(1:4,:));
fprintf(meshTxt,'</Nodes>\n<Elements type="hex8" mat="1" name="part1">\n');
elementsFormat = '<elem id=" %d "> %d, %d, %d, %d, %d, %d, %d, %d </elem>\n';
fprintf(meshTxt,elementsFormat,EleList2(1:9,:));
fprintf(meshTxt,'</Elements>\n</Geometry></febio_spec>');
fclose(meshTxt);

%%
figure
hold on
for i = 1:size(MeshBook,2)
    scatter3(MeshBook(1,i,:),MeshBook(2,i,:),i*ones(1,size(MeshBook,3)),'.','g');
    scatter3(MeshBook(1,1,:),MeshBook(2,1,:),i*ones(1,size(MeshBook,3)),'.','b');
end
hold off