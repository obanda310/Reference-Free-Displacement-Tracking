%% Radial Profiles
clear Radials Maxes
for i = 1:size(dirList,1)
    %try
    [MaxY, MaxX, RadY, RadX, radius] = RadialProfiles(dirList{i,1},prefix);
    radii(i,1) = radius*.1625;
    Maxes(1:size(MaxY,1)) = MaxY;
    Radials(1:size(RadY,1),1:size(RadY,2),i) = RadY;
    load([dirList{i,1} '\Shear Mat Files\ShearXYCutOff.mat'])
    cutOffs(i,1) = xyStd2;
    disp(num2str(i))
    %catch
    %   disp('Radial Profiles Failed')
    %  i
    %end
end



clear RadialsM RadialsM2
for i = 1:size(Radials,2)
    for j = 1:size(Radials,3)
        zeroI = find(Radials(4:end,i,j)<.3,1,'first');
        Radials(zeroI+3:end,i,j) = 0;
    end
end

%%
Radials(Radials == 0) = NaN;
Radials(Radials < mean(cutOffs)/1000) = NaN;
Radials(isnan(Radials)) = .3;
RadialsM = mean(Radials,3,'omitnan');
for i = 1:round(size(RadialsM,2)/2)
    RadialsM2(:,i) = mean([RadialsM(:,i),RadialsM(:,end-(i-1))],2,'omitnan');
end

%Pair up opposing slices (i.e. opposite angle relative to horizon)
for i = 1:round(size(Radials,2)/2)
    pairs(i,1) = i;
    pairs(i,2) = (size(Radials,2)+1)-i;
end

%Create colormap for associating profiles to pairs of slices
colormap = brewermap(size(pairs,1),'Dark2');



clear RadialsStd
for i = 1:size(RadialsM2,2)
    p1(:,:) = Radials(:,pairs(i,1),:);
    p2(:,:) = Radials(:,pairs(i,2),:);
    RadialsStd(:,i) = std([p1,p2],0,2,'omitnan') ;
end

RadialsEr = RadialsStd./sqrt(size(p1,2)*2);

cd(ListPath)
save('RadialData.mat')

%%
MeanProfile = figure;
hold on
Xs  = 0:.1625:(.1625*(size(RadialsM2,1)-1));
for i = 1:size(RadialsM2,2)
    plot(Xs,RadialsM2(:,i),'Color',colormap(i,:))
    %plot(Xs,RadialsM2(:,i)+RadialsStd(:,i),'Color',colormap(i,:))
    %plot(Xs,RadialsM2(:,i)-RadialsStd(:,i),'Color',colormap(i,:))
end
axis([0 25 0 6])
xt = 'Distance From Pattern (\mum)';
yt = 'Magnitude of Displacement (\mum)';
xl = xlabel(xt);
yl = ylabel(yt);
%tl = title(tt);
set(xl, 'fontweight','bold','fontsize',14);
set(yl,'fontweight','bold','fontsize',14);
title = ['\Mean ' prefix ' Plots'];
savefile = [ListPath title];
export_fig(MeanProfile,savefile,'-native');
cd(ListPath)
%%
%CHANGE FONT SIZES HERE
AxisFontSize = 28;
AxisTitleFontSize = 28;
LegendFontSize = 20;

AllProfiles = figure;
set(gcf,'position',[100,100,600,900])
hold on
for j = 1:size(Radials,2)
    [cidx, ~] = find(pairs == j, 1, 'first');
    for i = 1:size(Radials,3)
        plot(Xs,Radials(:,j,i),'Color',colormap(cidx,:))
    end
end
axis([0 25 0 7])
xt = {'Distance From'; 'Pattern (\mum)'};
yt = {'Magnitude of' ;'Displacement (\mum)'};
set(gca,'fontsize',AxisFontSize)
tt = 'Line-Profile Displacements';%input('enter the title','s');
le = 'Shear'; %input('enter the legend','s');
le2 = 'Normal';
le3 = 'Border';
xl = xlabel(xt);
yl = ylabel(yt);
%tl = title(tt);

set(xl, 'fontweight','bold','fontsize',AxisTitleFontSize);
set(yl,'fontweight','bold','fontsize',AxisTitleFontSize);
leg = legend([p1 p2 p3],le,le2,le3,'location','northwest');
leg.FontSize = LegendFontSize;
%set(tl,'fontweight','bold','fontsize',title_font_size)
title = ['\All ' prefix ' Plots'];
savefile = [ListPath title];
export_fig(AllProfiles,savefile,'-native');
cd(ListPath)

%% Plot Averages over individuals
colOptions{1,1} = 'white';
colOptions{2,1} = 'black';
colOptions{1,2} = 'black';
colOptions{2,2} = 'white';

for cI = 1:2 
        fcolor = colOptions{1,cI};
    bcolor = colOptions{2,cI};
    set(0,'defaultfigurecolor',bcolor)
%CHANGE FONT SIZES HERE
AxisFontSize = 40;
AxisTitleFontSize = 40;
LegendFontSize = 30;
Xs  = 0:.1625:(.1625*(size(RadialsM2,1)-1));
AllProfiles2 = figure;
set(gcf,'position',[100,100,1200,650])
hold on
% if strcmp(prefix ,'pill ') == 1
%     rectangle('Position',[0,0,14,7],'FaceColor',[.7 .7 .9020],'EdgeColor',[.7 .7 .9020])
% else
% rectangle('Position',[0,0,8,7],'FaceColor',[.65 .85 1],'EdgeColor',[.9 0 .9])
% end
    
    for j = 1:size(Radials,2)
        
        [cidx, ~] = find(pairs == j, 1, 'first');
        if strcmp(prefix ,'pill ') == 1
            for i = 1:size(Radials,3)
                p1 = plot(Xs,Radials(:,j,i),'Color',colormap(cidx,:),'LineWidth',1.5);
            end
            
        else
            for i = 1:size(Radials,3)
                p1 = plot(Xs,Radials(:,j,i),'Color',[.8 .8 .8],'LineWidth',1.5);
            end
        end
    end
    
    if strcmp(prefix ,'pill ') == 1
        for i = 1:size(RadialsM2,2)
            p2 = plot(Xs,RadialsM2(:,i),'Color',colormap(i,:),'LineWidth',9,'LineStyle','-.');
        end
    else
        
        p2 = plot(Xs,mean(RadialsM2,2,'omitnan')','Color',fcolor,'LineWidth',9,'LineStyle','-.');
              
    end
    
    %plot([0 20],[0.5 0.5],'LineWidth',3,'LineStyle','--','Color','black')
    
    
    
    if strcmp(prefix ,'pill ') == 1
        yticks([0:2:7])
        axis([0 20 0 7])
    else
        yticks([0:3])
        axis([0 20 0 3])
    end
    
    rectangle('Position',[0,0,25,.3],'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7])
    text(9,0.2,'Noise Range','FontSize',20)
    xt = {'Distance From Pattern (\mum)'};
    yt = {'Displacement (\mum)'}; %'Magnitude of' ;
    set(gca,'Color',bcolor)
    set(gca,'fontsize',AxisFontSize)
    set(gca,'fontsize',AxisFontSize,'XColor',fcolor,'YColor',fcolor,'LineWidth',2)
    tt = 'Line-Profile Displacements';%input('enter the title','s');
    le = 'Single Trace'; %input('enter the legend','s');
    le2 = 'Average Trace';
    le3 = 'Border';
    xl = xlabel(xt);
    yl = ylabel(yt);
    %tl = title(tt);
    
    set(xl, 'fontweight','bold','fontsize',AxisTitleFontSize,'color',fcolor);
    set(yl,'fontweight','bold','fontsize',AxisTitleFontSize,'color',fcolor);
    %leg = legend([p1 p2],le,le2,'location','Best');
    %leg.FontSize = LegendFontSize;
    %set(tl,'fontweight','bold','fontsize',title_font_size)
    title = ['\All ' prefix ' Plots 2 '  fcolor ' on ' bcolor];
    savefile = [ListPath title];
    export_fig(AllProfiles2,savefile,'-native');
    cd(ListPath)

end


%%
if strcmp(prefix ,'pill ') == 1
    RadialsM3 = RadialsM2(:,3)';
else
    RadialsM3 = mean(RadialsM2,2,'omitnan')';
end
RadialsM3(isnan(RadialsM3)) = [];
RadialsM4 = RadialsM3-0.5;
Xs2 = [0:0.1625:(0.1625*(size(RadialsM4,2)-1))];
RadM4 = fit(Xs2',RadialsM4','cubicinterp');
figure
plot(RadM4)
if strcmp(prefix ,'pill ') == 1
    fzero(RadM4,14)
else
    fzero(RadM4,7)
end

%%
Xs3 = [0:0.1625:(0.1625*(size(Radials,1)-1))];
if strcmp(prefix ,'pill ') == 1
    for i = 1:size(Radials,3)
        for j = 1:size(Radials,2)
            fun1 = fit(Xs3',Radials(:,j,i)-.5,'cubicinterp');
            zeros(i,j) = fzero(fun1,10);
        end
    end
    
    for i = 1:size(pairs,1)
        RadialsZeros(:,i) = cat(1,zeros(:,pairs(i,1)),zeros(:,pairs(i,2)));
    end
    
    RadialsZS(1,1:3) = mean(RadialsZeros);
    RadialsZS(2,1:3) = std(RadialsZeros);
    save('RadialZerosPill.mat','RadialsZS','zeros')
else
    figure
    hold on
    
    RadialsC = Radials - .5;
    for i = 1:size(RadialsC,3)
        
        
        
        
        
        for j = 1:size(RadialsC,2)
            fun1 = fit(Xs3',RadialsC(:,j,i),'cubicinterp');
            %plot(fun1)
            %xlim([0,10])
            %ylim([-1,10])
            try
                idx = find(RadialsC(:,j,i)>0,1,'last');
                zeros(i,j) = fzero(fun1,[idx]);
                
            catch
                zeros(i,j) = NaN;
            end
            %            plot([zeros(i,j) zeros(i,j)],[0 1])
            
        end
    end
    zeros(zeros>20) = 0;
    zeros(isnan(zeros)) = 0;
    clear RadialsZS
    RadialsZS(1,1) = mean(zeros(:));
    RadialsZS(2,1) = std(zeros(:));
    save('RadialZerosCircle.mat','RadialsZS','zeros')
    
end
