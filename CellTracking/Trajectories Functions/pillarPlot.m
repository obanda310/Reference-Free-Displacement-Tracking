function trajOfInterest = pillarPlot(book1,book2,book3,cm3,imageBlack,totalNumFrames,interpBook,totalAverageInterp)
cm3(cm3 == 0) = NaN;

        plottedProfiles = figure('units','pixels','name','Select Outputs','position',[100,100,1000,750],...
    'toolbar','none','menu','none');
% Starting point for pillar fitting
trajOfInterest = 1;
handles = guidata(gcf);
handles.checkbox = uicontrol('style','checkbox','units','normalized',...
    'position',[.10,.150,.100,.03],'string','Attempt Fit','Value',0);

trajInt()



handles.trajOfInt = uicontrol( ...
    'Style','Edit', ...
    'String','1', ...
    'Units','normalized', ...
    'Position',[.100,.250,.100,.030], ...
    'Callback',@trajIntListener);

 function trajInt(~,~)
     
        interpXq = linspace(1,size(book1,2),size(book1,2)*3);
        
        subplot(2,3,1:3)
        plot(linspace(1,totalNumFrames,totalNumFrames),book1(7,:,trajOfInterest),'.','MarkerSize',5,'Color',[0 0 0])
        hold on
        plot(interpXq,interpBook(trajOfInterest,:),'.','MarkerSize',5,'Color',[0 0 0])
        
        %create a smoothed curve for finding maxima and minima
        trajOfIntInt = book1(6,:,trajOfInterest);
        trajOfIntInt(trajOfIntInt == 0) = NaN;
        trajOfIntIntMean = trajOfIntInt;
        for i = 3:size(trajOfIntInt,2)-2
            trajOfIntIntMean(1,i) = mean(trajOfIntInt(1,i-2:i+2));
        end
        
        plot(linspace(1,totalNumFrames,totalNumFrames),trajOfIntIntMean,'.','MarkerSize',20,'Color',[0 0 0])
        
        findpeaks((book1(6,:,trajOfInterest)))
        findpeaks((interpBook(trajOfInterest,:)))
        findpeaks(totalAverageInterp(1,:))
        
        currentAverages = zeros(totalNumFrames,1);
        
        for i = 1:totalNumFrames
            errorbar(i,mean(cm3(trajOfInterest,1:50,i),'omitnan'),std(cm3(trajOfInterest,1:50,i),'omitnan'),'.','MarkerSize',5,'Color','r');
            currentAverages(i,1) = mean(cm3(trajOfInterest,1:50,i),'omitnan');
        end
        currentAverages(isnan(currentAverages))=0;
        interpCurrentAverages = interp1(linspace(1,size(book1,2),size(book1,2)),conv(currentAverages(:,1),gausswin(6),'same'),interpXq,'spline');
        findpeaks(interpCurrentAverages(1,:))
        
        
        %fit the traj of interest with pillarFit.m
        vals = p_call();
        if ismember(1,vals) == 1
        attemptPillarFit(currentAverages)
        end
        hold off
        
        
        %plot the locations of traj of interest and surrounding undeformed
        %trajs
        subplot(2,3,5)
        imshow(imageBlack,[])
        hold on
        for i = 1:size(book3,2)
            scatter(book2(book3(trajOfInterest,i),1),book2(book3(trajOfInterest,i),2),'.','g')
        end
        scatter(book2(trajOfInterest,1),book2(trajOfInterest,2),'o','b')
        
        
        %plot magnitude of displacement for the traj of interest
        subplot(2,3,6)
        
        plot(linspace(1,totalNumFrames,totalNumFrames),book1(5,:,trajOfInterest),'.','MarkerSize',5,'Color',[0 0 0])
        axis([0,totalNumFrames,0,max(max(max(book1(5,:,:))))])

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