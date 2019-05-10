function [Surface2,SurfaceAll,Zeros] = findSurface2(shear,image,raw)

RawStack = permute(image.RawStack,[2,1,3]);
originalTest = image.RawStack(:,:,end);

oMean = median(originalTest(:));
oStd = std(originalTest(originalTest(:)<prctile(originalTest(:),80)));

intCutoff = oMean+2*oStd;

%% Pull intensity values from the raw image stack
% Should be much more accurate than using the processed stack for surface
% finding.

clear Int2
for i = 1:size(shear.Int,2)
    for j = 1:size(RawStack,3)
        if j>size(shear.Int,1)
            X = round(shear.rawX(shear.lastFrame(i),i)/raw.dataKey(9,1));
            Y = round(shear.rawY(shear.lastFrame(i),i)/raw.dataKey(9,1));
            
        elseif (shear.rawX(j,i))>0
            X = round(shear.rawX(j,i)/raw.dataKey(9,1));
            Y = round(shear.rawY(j,i)/raw.dataKey(9,1));
        else
            %disp('Zeros!')
            X = round(shear.rawX(shear.lastFrame(i),i)/raw.dataKey(9,1));
            Y = round(shear.rawY(shear.lastFrame(i),i)/raw.dataKey(9,1));
        end
        if X>1 && X<size(RawStack,2) && Y>1 && Y<size(RawStack,1)
            tempint = RawStack(Y-1:Y+1,X-1:X+1,j);
            Int2(j,i) = mean(tempint(:))-intCutoff;
            %         else
            %             Int2(j,i) = RawStack(Y,X)-intCutoff;
        end
    end
end
%%
Int3 = Int2;
Int3(Int3<0) = 0;
for i = 1:size(Int3,2)
    IntInd(i) = find(Int3(:,i),1,'last');
    if  IntInd(i)>3         
        if(Int3(IntInd(i),i) + Int3(IntInd(i)-1,i) + Int3(IntInd(i)-2,i)) == Int3(IntInd(i),i) %if last signal is isolated (probably noise)
            cont = 1;
            while cont == 1
            Int3(IntInd(i),i) = 0; %then get rid of it
            IntInd(i) = find(Int3(:,i),1,'last'); %update last signal
            if (Int3(IntInd(i),i) + Int3(IntInd(i)-1,i) + Int3(IntInd(i)-2,i)) == Int3(IntInd(i),i) %Continue checking until real end is found
                Int3(IntInd(i),i) = 0; %then get rid of it
            IntInd(i) = find(Int3(:,i),1,'last'); %update last signal
            else
                cont = 0;
            end
            end       
        end        
    end
end

topSurface = [0,0,0];
surfaceFilt = image.ADil&image.Borders;
for i = 1:size(IntInd,2)
    %if it is under the cell
    
    if surfaceFilt(round((shear.lastY(i))/raw.dataKey(9,1)),round((shear.lastX(i))/raw.dataKey(9,1)))~=0
        topSurface = cat(1,topSurface,[shear.lastX(i) shear.lastY(i) IntInd(i)*raw.dataKey(10,1)]);
    end
end
%shift cells up 1 to get rid of initial zero
topSurface(1,:) = [];
Surface2 = fit([topSurface(:,1),topSurface(:,2)],topSurface(:,3),'poly11');

Zeros(:,1) = shear.lastX;
Zeros(:,2) = shear.lastY;
Zeros(:,3) = IntInd;

SurfaceAll = fit([Zeros(:,1),Zeros(:,2)],Zeros(:,3),'poly11');
%% Show surface features and normal features
figure
scatter3(shear.lastX,shear.lastY,IntInd)
xlim([0 max(shear.rawX(:))])
ylim([0 max(shear.rawY(:))])
zlim([0 size(RawStack,3)*raw.dataKey(10,1)])
hold on
%scatter3(r.r(:,1),r.r(:,2),r.r(:,3))
plot(Surface2)
plot(SurfaceAll)

%%
profiles = figure;
map = brewermap(size(Int2,2),'spectral');
hold on
for i = 1:size(Int2,2)
    plot3(1:1:size(RawStack,3),Int2(:,i),i*ones(1,size(RawStack,3)),'color',map(i,1:3))
end
patch([0 size(RawStack,3) size(RawStack,3) 0], [0 0 0 0],[0 0 size(Int2,2) size(Int2,2)],'black')
view([180 -45])
savefile = 'PillarIntProfiles.tif';
export_fig(profiles,savefile,'-native');
% 
% %% Fit intensity profiles and find where they become noise
% %(i.e. intersect the surface)
% clear failcount
% failcount = 0;
% Int2 = double(Int2);
% progressbar('Approximating Substrate Surface')
% for i = 1:size(shear.rawZ,2)
%     lastZ(i) = round(shear.rawZ(shear.lastFrame(1,i),i)/raw.dataKey(10,1));
%     fitStart(i) = lastZ(i) - 2;
% end
% % meanZ = floor(mean((lastZ)))
% % stdZ = 3*ceil(std(lastZ))
% %fitEnd = meanZ - stdZ;
% fitEnd1 = zeros(size(fitStart));
% for i = 1:size(Int2,2)
%     progressbar(i/size(Int2,2))        
%     
%     %Create an endpoint for zero search window
%     try
%         fitEnd1(i) = lastZ(i) + find(Int2(lastZ(i):end,i)<0,1,'first')+2;
%         tempInt = Int2(fitStart(i):fitEnd1(i),i);
%         Xs = (fitStart(i):fitEnd1(i))';
%         %tempInt(tempInt==0)=NaN;
%         fun1 = fit(Xs,tempInt,'smoothingspline');
%     catch
%         fitEnd1(i) = size(RawStack,3);
%         tempInt = Int2(fitStart(i):fitEnd1(i),i);
%         Xs = (fitStart(i):fitEnd1(i))';
%         %tempInt(tempInt==0)=NaN;
%         fun1 = fit(Xs,tempInt,'smoothingspline');
%     end
%     
%     
% %     found = 0;
% %     count = 0;
% %     nozero = 0;
% %     maxF = size(RawStack,3);
% %     
% %     %Determine if a zero exists on interval
% %     while found == 0        
% %         funcheck = feval(fun1,fitStart(i)+count);        
% %         if funcheck <0
% %             fitEnd1(i) = fitStart(i)+count;
% %             found = 1;
% %         end        
% %         if fitStart(i)+count > maxF
% %             found = 1;
% %             nozero = 1;
% %         end
% %         count = count+1;       
% %     end
% %     
% %     
% %     if nozero == 1
% %         Zeros(i,3) = NaN;
% %     else
% 
%         %Find Zero within search window
%         try
%             Zeros(i,3) = fzero(fun1,[fitStart(i) fitEnd1(i)]);
%         catch
%             failcount = failcount+1;
%             %             [fitStart fitEnd]
%             %             feval(fun1, fitStart)
%             %             feval(fun1, fitEnd)
%         end
% %     end
%     
%     Zeros(i,1) = shear.rawX(shear.lastFrame(i),i);
%     Zeros(i,2) = shear.rawY(shear.lastFrame(i),i);    
%     
%     %plot(Xs,tempInt)    
%     %scatter(zeros, feval(fun1,zeros))
%     %fun1 = fit(Xs,tempInt,'cubicinterp');    
% end
% %%
% for i = 1:size(fitStart,2)
%     Int3(1,i) = Int2(fitStart(i),i);
% end
% 
% scatter3(Zeros(:,3),zeros(size(Zeros(:,3))),[1:(size(Zeros(:,3),1))]','filled','MarkerFaceColor','red');
% scatter3(fitStart,Int3,[1:(size(Zeros(:,3),1))]','filled','MarkerFaceColor','blue');
% scatter3(fitEnd1,zeros(size(Zeros(:,3))),[1:(size(Zeros(:,3),1))]','filled','MarkerFaceColor','green');
% Zeros(:,3) = Zeros(:,3)*raw.dataKey(10,1);
% 
% %% Show surface features and normal features
% % figure
% % scatter3(zeros(:,1),zeros(:,2),zeros(:,3))
% % xlim([0 max(shear.rawX(:))])
% % ylim([0 max(shear.rawY(:))])
% % zlim([0 size(RawStack,3)*raw.dataKey(10,1)])
% % hold on
% % scatter3(r.r(:,1),r.r(:,2),r.r(:,3))
% % try
% % scatter3(emptyZ(:,1), emptyZ(:,2), emptyZ(:,3),'b')
% % end
% %% Fit planar surface to non-deformed dots
% % The 'shear' script already did basically this, just need to repeat using
% % new points.
% clear zerosnNaN
% topSurface = [0,0,0];
% surfaceFilt = image.ADil&image.Borders;
% for i = 1:size(Zeros,1)
%     %if it is under the cell
%     if surfaceFilt(round((Zeros(i,2))/raw.dataKey(9,1)),round((Zeros(i,1))/raw.dataKey(9,1)))~=0
%         topSurface = cat(1,topSurface,Zeros(i,1:3));
%     end
% end
% %shift cells up 1 to get rid of initial zero
% topSurface(1,:) = [];
% 
% %remove outliers (typically incomplete pillars at the edges)
% ub = mean(topSurface(:,3)) + 2*std(topSurface(:,3));
% lb = mean(topSurface(:,3)) - 2*std(topSurface(:,3));
% topSurface([find(topSurface(:,3)>ub | topSurface(:,3)<lb)],:) = [];
% topSurface(isnan(topSurface(:,3)),:) = [];
% zerosNaN = find(isnan(Zeros(:,3)));
% zerosnNaN(1:size(Zeros,1),1) = 1:size(Zeros,1);
% zerosnNaN(zerosNaN,:) = [];
% 
% 
% %Surface2 = fit([topSurface(:,1),topSurface(:,2)],topSurface(:,3),'lowess','Span',0.1);
% Surface3 = fit([topSurface(:,1),topSurface(:,2)],topSurface(:,3),'poly11');
% SurfaceAll = fit([Zeros(zerosnNaN,1),Zeros(zerosnNaN,2)],Zeros(zerosnNaN,3),'lowess','Span',0.005);
% 
% 
% %% Show surface features and normal features
% figure
% scatter3(Zeros(:,1),Zeros(:,2),Zeros(:,3))
% xlim([0 max(shear.rawX(:))])
% ylim([0 max(shear.rawY(:))])
% zlim([0 size(RawStack,3)*raw.dataKey(10,1)])
% hold on
% %scatter3(r.r(:,1),r.r(:,2),r.r(:,3))
% plot(Surface2)
% plot(SurfaceAll)

%%

toc