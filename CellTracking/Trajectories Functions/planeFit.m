function planeFit(book1,book2,totalNumFrames,filePath)
 close all
    numTraj = size(book1,3);
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
    
    folderName = ('Planefit');
    mkdir(filePath,folderName)
    filePathPlane = [filePath '\Planefit'];
    
    savefile = [filePathPlane '\Intensity Profile Linear Fits.tif'];
    export_fig(intProfilesFits,savefile);
    hold off
    
    intProfilesRaw = figure;
    image(intProfiles,'cdatamapping','scaled'); colormap(gray);
    savefile = [filePathPlane '\Intensity Profiles 1 Raw.tif'];
    export_fig(intProfilesRaw,savefile);
    
    intProfilesPSfig = figure;
    image(intProfilesPS,'cdatamapping','scaled'); colormap(gray);
    savefile = [filePathPlane '\Intensity Profiles 2 Phase Sym.tif'];
    export_fig(intProfilesPSfig,savefile);
    
    intProfilesMaxfig = figure;
    image(intProfilesMax,'cdatamapping','scaled'); colormap(gray);
    savefile = [filePathPlane '\Intensity Profiles 3 Max.tif'];
    export_fig(intProfilesMaxfig,savefile);
    
    intProfilesClosedfig1 = figure;
    image(intProfilesClosed1,'cdatamapping','scaled'); colormap(gray);
    savefile = [filePathPlane '\Intensity Profiles 4 Closed.tif'];
    export_fig(intProfilesClosedfig1,savefile);
    
    intProfilesCleanedfig = figure;
    image(intProfilesCleaned,'cdatamapping','scaled'); colormap(gray);
    savefile = [filePathPlane '\Intensity Profiles 5 Cleaned.tif'];
    export_fig(intProfilesCleanedfig,savefile);
    
    intProfilesClosedfig2 = figure;
    image(intProfilesClosed2,'cdatamapping','scaled'); colormap(gray);
    savefile = [filePathPlane '\Intensity Profiles 6 Closed.tif'];
    export_fig(intProfilesClosedfig2,savefile);
    
    intProfilesLabeledfig = figure;
    image(intProfilesLabeled,'cdatamapping','scaled'); colormap(gray);
    savefile = [filePathPlane '\Intensity Profiles 7 Labeled.tif'];
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