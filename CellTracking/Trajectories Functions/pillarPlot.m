function trajOfInterest = pillarPlot(book1,book2,book3,cm3,roiStack,totalNumFrames,interpBook,totalAverageInterp)
cm3(cm3 == 0) = NaN;

plottedProfiles = figure('units','pixels','name','Visually Analyzing Z-Axis Deformation and Compression','position',[100,100,1000,750],...
    'toolbar','none','menu','none');

% Starting point for pillar fitting
trajOfInterest = 1;
handles = guidata(gcf);
handles.checkbox = uicontrol('style','checkbox','units','normalized',...
    'position',[.10,.150,.100,.03],'string','Attempt Fit','Value',0);

handles.origPerTextBox = uicontrol('style','text','Units','normalized','position',[.100,.210,.150,.030]);
        origPeriod = 0;
        origStr = num2str(origPeriod);
        set(handles.origPerTextBox,'String',['Current Pillar Period: ' origStr ''])
        
        handles.avgPerTextBox = uicontrol('style','text','Units','normalized','position',[.100,.180,.150,.030]);
        avgPeriod = 0;
        avgStr = num2str(avgPeriod);
        set(handles.avgPerTextBox,'String',['Average Pillar Period: ' avgStr ''])

trajInt()



handles.trajOfInt = uicontrol( ...
    'Style','Edit', ...
    'String','1', ...
    'Units','normalized', ...
    'Position',[.100,.250,.100,.030], ...
    'Callback',@trajIntListener);

    function trajInt(~,~)
        
        %x values for plots
        frameXs = linspace(1,totalNumFrames,totalNumFrames);
        interpXq = linspace(1,size(book1,2),size(book1,2)*3);
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
%%%%%%%%%%%%%%%%%Top Plot%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

        subplot(4,3,1:3)
        hold off
        %Plots Original Intensity Values
        plot(frameXs,book1(6,:,trajOfInterest),'Color',[0 0 .7]) 
        axis([0,totalNumFrames,0,max(max(max(book1(6,:,:))))])
        hold on
        currentAverages = zeros(totalNumFrames,1);
        for i = 1:totalNumFrames
            %plots stdev from mean intensity
            errorbar(i,mean(cm3(trajOfInterest,1:50,i),'omitnan'),std(cm3(trajOfInterest,1:50,i),'omitnan'),'.','MarkerSize',5,'Color',[0 .7 0]); 
            currentAverages(i,1) = mean(cm3(trajOfInterest,1:50,i),'omitnan');
        end
        currentAverages(isnan(currentAverages))=0;
        
        title('Compare Current Pillar to Nearby Average Pillar')
        ylabel('Intensity')
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%%%%%%%%%%%%%%%%Middle Plot%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        hold off
        subplot(4,3,4:6)       
        
        %plot original points with gaussian applied
        plot(frameXs,book1(7,:,trajOfInterest),'o','MarkerSize',5,'Color',[0 0 .7])
        hold on
        [gaussPeaks,gaussLocs] = findpeaks((book1(7,:,trajOfInterest)));
        clear realPeaks
        realPeaks = find(gaussPeaks>max(gaussPeaks)/4);
        plot(frameXs(gaussLocs(realPeaks)),gaussPeaks(realPeaks),'x','MarkerSize',15,'Color',[.5 0 0]) 
        
        %plot interpolated intensity values of current pillar
        plot(interpXq,interpBook(trajOfInterest,:),'.','MarkerSize',5,'Color',[0 0 0])
        [interpPeaks,interpLocs] = findpeaks(interpBook(trajOfInterest,:));
        clear realPeaks
        realPeaks = find(interpPeaks>max(interpPeaks)/4);
        plot(interpXq(interpLocs(realPeaks)),interpPeaks(realPeaks),'x','MarkerSize',15,'Color',[0 0 .7]) 
        axis([0,totalNumFrames,0,max(max(interpBook))])
        origPeriod = (mean(interpLocs(realPeaks(3:end-1))-interpLocs(realPeaks(2:end-2)))/3)*.4;        
        handles.origPerTextBox = uicontrol('style','text','Units','normalized','position',[.100,.210,.15,.030]);
        origStr = num2str(origPeriod);
        set(handles.origPerTextBox,'String',['Current Pillar Period: ' origStr ''])
        
        %plot original data with peaks visible
        plot(frameXs,book1(6,:,trajOfInterest),'Color',[0 0 .7]) 
        [origPeaks,origLocs] = findpeaks((book1(6,:,trajOfInterest)));
        plot(frameXs(origLocs),origPeaks,'x','MarkerSize',15,'Color',[0 0 .7])
        title('Compare Raw Pillar to Processed Pillar')
        ylabel('Intensity')
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%Bottom Plot%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        hold off
        subplot(4,3,7:9) %Hold is off so plot is cleared on refresh
        
        plot(interpXq,interpBook(trajOfInterest,:),'Color',[0 0 .7])
        hold on
        [origInterpPeaks,origInterpLocs] = findpeaks((interpBook(trajOfInterest,:)));
        clear realPeaks
        realPeaks = find(origInterpPeaks>max(origInterpPeaks)/4);
        plot(interpXq(origInterpLocs(realPeaks)),origInterpPeaks(realPeaks),'o','MarkerSize',10,'Color',[0 0 .7])
        
        
        
%         plot(interpXq,totalAverageInterp(1,:))
%         [taiPeaks,taiLocs] = findpeaks(totalAverageInterp(1,:));
%         plot(interpXq(taiLocs),taiPeaks,'x','MarkerSize',15,'Color',[.5 0 0])
        
        interpCurrentAverages = interp1(linspace(1,size(book1,2),size(book1,2)),conv(currentAverages(:,1),gausswin(6),'same'),interpXq,'spline');
        plot(interpXq,interpCurrentAverages(1,:),'Color',[0 .8 0])
        [icaPeaks,icaLocs] = findpeaks(interpCurrentAverages(1,:));
        clear realPeaks
        realPeaks = find(icaPeaks>max(icaPeaks)/4);
        plot(interpXq(icaLocs(realPeaks)),icaPeaks(realPeaks),'o','MarkerSize',10,'Color',[0 .8 0])
        meanPeriod = (mean(interpLocs(realPeaks(2:end-1))-interpLocs(realPeaks(1:end-2)))/3)*.4;
        handles.avgPerTextBox = uicontrol('style','text','Units','normalized','position',[.100,.180,.150,.030]);
        avgPeriod = meanPeriod;
        avgStr = num2str(avgPeriod);
        set(handles.avgPerTextBox,'String',['Average Pillar Period: ' avgStr ''])
        hold off
        
        title('Compare Processed Pillar to Processed Average')
        xlabel('Frame')
        ylabel('Intensity')
        axis([0,totalNumFrames,0,max(max(interpBook))])
        
        %fit the traj of interest with pillarFit.m
        vals = p_call();
        if ismember(1,vals) == 1
            attemptPillarFit(currentAverages)
        end

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%Bottom Right Plot%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

        %plot magnitude of displacement for the traj of interest
        subplot(4,3,12)
        plot(linspace(1,totalNumFrames,totalNumFrames),book1(5,:,trajOfInterest),'.','MarkerSize',5,'Color',[0 0 0])
        axis([0,totalNumFrames,0,max(max(max(book1(5,:,:))))])
        title('Pillar Displacement Through Z')
        xlabel('Frame')
        ylabel('Displacement[um]')
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%Image Stack Plot%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

        %plot the locations of traj of interest and surrounding undeformed
        %trajs
        
        ShowStack(roiStack,438)
        hold on
        for i = 1:size(book3,2)
            scatter(book2(book3(trajOfInterest,i),1),book2(book3(trajOfInterest,i),2),'.','g')
        end
        scatter(book2(trajOfInterest,1),book2(trajOfInterest,2),'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[1 1 1],'LineWidth',2)
        hold off
                
        
    end
    
 

    function trajIntListener(~,~)
        
        trajOfInterest = str2double(get(handles.trajOfInt,'String'));
        trajInt()
    end

    function attemptPillarFit(currentAverages,~)
        hold on
        xVals = linspace(1,totalNumFrames,totalNumFrames);
        fitWeights = zeros(totalNumFrames,1);
        fitWeights(:,1) = xVals(1,:)+1000;
        s = fitoptions('Method','NonlinearLeastSquares',...
            'Lower',[4000,25000,0,0,0,0,0],...
            'Upper',[6000,35000,15,.1,1,3.14,15],...
            'StartPoint',[5500,29000,11.5,.015,0.1731,2.2,6.712]);
        f = fittype('pillarFit2(x,S,iS,iP,P,Pg,Mo,Mf)','options',s);
        [fitCoef,~] = fit(xVals',currentAverages,f,'weight',fitWeights);
        fitCoef2 = coeffvalues(fitCoef);
        yVals = pillarFit(xVals,fitCoef2(1,5),fitCoef2(1,7),fitCoef2(1,6),fitCoef2(1,3),fitCoef2(1,4),fitCoef2(1,2),fitCoef2(1,1));
        plot(xVals,yVals,'MarkerSize',20,'Color','bl')
        findpeaks(yVals(1,:))
        
        hold off
    end

    function vals = p_call(varargin)
        if isempty(get(handles.checkbox,'Value'))
            vals = 0;
        else
            vals = get(handles.checkbox,'Value');
        end
    end
end