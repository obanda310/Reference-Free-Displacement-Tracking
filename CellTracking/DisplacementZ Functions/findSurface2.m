function [Surface2,SurfaceAll,Zeros] = findSurface2(shear,image,raw)

image.RawStack = permute(image.RawStack,[2,1,3]);
originalTest = image.RawStack(:,:,35);
oMean = mean(originalTest(:));
oStd = std(originalTest(:));
intCutoff = oMean+2*oStd;

%% Pull intensity values from the raw image stack
% Should be much more accurate than using the processed stack for surface
% finding.

clear Int2
for i = 1:size(shear.Int,2)
    for j = 1:size(image.RawStack,3)
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
        if X>1 && X<size(image.RawStack,2) && Y>1 && Y<size(image.RawStack,1)
        tempint = image.RawStack(Y-1:Y+1,X-1:X+1,j);
        Int2(j,i) = mean(tempint(:))-intCutoff;
%         else
%             Int2(j,i) = image.RawStack(Y,X)-intCutoff;
        end
    end
end
%%
figure
hold on
for i = 1:size(Int2,2)
    plot(Int2(:,i))
end
plot([0 size(image.RawStack,3)], [0 0])
    

%% Fit intensity profiles and find where they become noise
%(i.e. intersect the surface)
Int2 = double(Int2);
progressbar('Approximating Substrate Surface')
for i = 1:size(Int2,2)
    progressbar(i/size(Int2,2))
    fitStart = shear.lastFrame(i)-5;
    fitEnd1 = shear.lastFrame(i) + find(Int2(shear.lastFrame(i):end,i)<0,1,'first')+2;
    tempInt = Int2(fitStart:fitEnd1,i);
    Xs = (fitStart:fitEnd1)';
    tempInt(tempInt==0)=NaN;
    fun1 = fit(Xs,tempInt,'smoothingspline');
    
    found = 0;
    count = 0;
    nozero = 0;
    maxF = size(image.RawStack,3);
    while found == 0 
        
        funcheck = feval(fun1,fitStart+count);
        
        if funcheck <0
            fitEnd = fitStart+count;
            found = 1;
        end
        
        if fitStart+count > maxF
            found = 1;
            nozero = 1;
        end
        count = count+1;
        
    end
    
    if nozero == 1
    Zeros(i,3) = NaN;    
    else
    Zeros(i,3) = fzero(fun1,[fitStart fitEnd]);
    end
    
    Zeros(i,1) = shear.rawX(shear.lastFrame(i),i);
    Zeros(i,2) = shear.rawY(shear.lastFrame(i),i);
    
   
    %plot(Xs,tempInt)
    
    %scatter(zeros, feval(fun1,zeros))
    %fun1 = fit(Xs,tempInt,'cubicinterp');
    
end
Zeros(:,3) = Zeros(:,3)*raw.dataKey(10,1);
%% Show surface features and normal features
% figure
% scatter3(zeros(:,1),zeros(:,2),zeros(:,3))
% xlim([0 max(shear.rawX(:))])
% ylim([0 max(shear.rawY(:))])
% zlim([0 size(image.RawStack,3)*raw.dataKey(10,1)])
% hold on
% scatter3(r.r(:,1),r.r(:,2),r.r(:,3))
% try
% scatter3(emptyZ(:,1), emptyZ(:,2), emptyZ(:,3),'b')
% end
%% Fit planar surface to non-deformed dots
% The 'shear' script already did basically this, just need to repeat using
% new points.
clear zerosnNaN
topSurface = [0,0,0];
surfaceFilt = image.ADil&image.Borders;
for i = 1:size(Zeros,1)
    %if it is under the cell
    if surfaceFilt(round((Zeros(i,2))/raw.dataKey(9,1)),round((Zeros(i,1))/raw.dataKey(9,1)))~=0
            topSurface = cat(1,topSurface,Zeros(i,1:3));        
    end
end
%shift cells up 1 to get rid of initial zero
topSurface(1,:) = [];

%remove outliers (typically incomplete pillars at the edges)
ub = mean(topSurface(:,3)) + 2*std(topSurface(:,3));
lb = mean(topSurface(:,3)) - 2*std(topSurface(:,3));
topSurface([find(topSurface(:,3)>ub | topSurface(:,3)<lb)],:) = [];
topSurface(isnan(topSurface(:,3)),:) = [];
zerosNaN = find(isnan(Zeros(:,3)));
zerosnNaN(1:size(Zeros,1),1) = 1:size(Zeros,1);
zerosnNaN(zerosNaN,:) = [];


%Surface2 = fit([topSurface(:,1),topSurface(:,2)],topSurface(:,3),'lowess','Span',0.1);
Surface2 = fit([topSurface(:,1),topSurface(:,2)],topSurface(:,3),'poly11');
SurfaceAll = fit([Zeros(zerosnNaN,1),Zeros(zerosnNaN,2)],Zeros(zerosnNaN,3),'lowess','Span',0.005);


%% Show surface features and normal features
% figure
% scatter3(Zeros(:,1),Zeros(:,2),Zeros(:,3))
% xlim([0 max(shear.rawX(:))])
% ylim([0 max(shear.rawY(:))])
% zlim([0 size(image.RawStack,3)*raw.dataKey(10,1)])
% hold on
% scatter3(r.r(:,1),r.r(:,2),r.r(:,3))
% %plot(Surface2)
% plot(SurfaceAll)

%%

toc