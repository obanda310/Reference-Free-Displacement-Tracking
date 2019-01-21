%% Automated Script
clear all
close all
[ListName,ListPath] = uigetfile;
fullListName = strcat(ListPath,ListName);
load(strcat(ListPath,ListName),'dirList')
parts = strsplit(ListName, 'List');
prefix = parts{1};
set(0,'defaultfigurecolor',[1 1 1])
%% For updating paths (when files are moved)
%set a as new path root
%set b as old path root
ak = strfind(fullListName,ListName);

a = [fullListName(1:ak-1) ''];

Schar = char(dirList{:});
all_rows_same = all(diff(Schar, [], 1) == 0, 1);
common_cols = find(~all_rows_same, 1, 'first');
if isempty(common_cols)
    b = '?'
else
    b = dirList{1}(1:common_cols-1);
end

for i = 1:size(dirList,1)
    
    dirList{i,1} = strrep(dirList{i,1},b,a);
end
%%
save(fullListName,'dirList')

%% For clearing old files
for j = 1:2
    for i = 1:size(dirList,1)
        cleanOld(dirList{i,1})
        disp(num2str(i))
    end
end

%% Create Variables for Trace Profiles
for i = 1:size(dirList,1)
    disp(dirList{i,1})
    cd(dirList{i,1});
    load('3Ddata.mat')
    filePath = cd;
    folderName = 'Profile Data';
    mkdir(filePath,folderName)
    save('Profile Data\vqZ.mat','vqN','vq3','image','HeatMapN','HeatMap3')
end


%% For Running Profile Script
clear caught goodData
goodData = [1;2;3;5;6;14;17;23;31;34;43;44;47;49;50;48;11;12;19;25;28;36];
for i = 1:size(goodData,1)
    try
        [normXY,normZ,normAxis,cell_boundary] = profFunc(dirList{goodData(i,1),1});
        profBook(1,1:size(normAxis,2),i) = normAxis;
        profBook(2,1:size(normXY,2),i) = normXY;
        profBook(3,1:size(normZ,2),i) = normZ;
        cb(i) = cell_boundary;
        save('Profile Data\ProfileOutputData.mat','normXY','normZ','normAxis')
    catch
        caught(i) = 1;
        i
        dirList{goodData(i,:),1}
    end
end
cbf = mean(cb);
save(strcat(ListPath,prefix,'Profiles.mat'),'profBook','cb','cbf')

%%
prefix = '\Gels';
cb(cb<25 | cb>75)=NaN;
cbf2 = mean(cb,'omitnan')/100;
aScale = 1/cbf2;
clear profBookCat
profBookCat = profBook;
clear keep keep2 keep3 map profBook2 profBook3
for i = 1:size(profBookCat,3)
    %max(profBookCat(2,:,i))
    if max(profBookCat(2,:,i))>1
        keep(1,i) = 1;
        keep(2,i) = max(profBookCat(2,:,i));
        profBook2(1:size(profBookCat(:,:,i),1),1:size(profBookCat(:,:,i),2),i) = profBookCat(:,:,i);
    end
end
profBook2(:,:,~keep(1,:)) = [];
keep2 = keep(2,:);
keep2(:,~keep(1,:)) = [];
[keep3,sortIdx] = sort(keep2);

colOptions{1,1} = 'white';
colOptions{2,1} = 'black';
colOptions{1,2} = 'black';
colOptions{2,2} = 'white';
map = brewermap(size(profBook2,3),'*spectral');

for i = 1:size(colOptions,2)
    fcolor = colOptions{1,i};
    bcolor = colOptions{2,i};
    
    set(0,'defaultfigurecolor',bcolor)
    ProfileOverlays = figure;
    hold on
    for i = 1:size(profBook2,3)
        plot(profBook2(1,:,sortIdx(i))*aScale,profBook2(2,:,sortIdx(i)),'color',[map(i,1:3)])
        plot(profBook2(1,:,sortIdx(i))*aScale,profBook2(3,:,sortIdx(i)),'color',[map(i,1:3)],'linestyle','-.')
    end
    p1= plot(profBook2(1,:,sortIdx(i))*aScale,mean(profBook2(2,:,sortIdx(:)),3,'omitnan'),'color',fcolor,'linewidth',3);
    p2=plot(profBook2(1,:,sortIdx(i))*aScale,mean(profBook2(3,:,sortIdx(:)),3,'omitnan'),'color',fcolor,'linestyle','-.','linewidth',3);
    p3=plot([cbf2 cbf2]*aScale,[min(min(min(profBook2))),max(max(max(profBook2)))],'color',[.5 .5 .5],'linestyle','--','linewidth',1);
    p4=plot([.4 .4],[min(min(min(profBook2))),2.2],'color',[1,.5,.5],'linestyle','--','linewidth',1);
    p5=plot([1.6 1.6],[min(min(min(profBook2))),max(max(max(profBook2)))],'color',[1,.5,.5],'linestyle','--','linewidth',1);
    set(gca,'Color',bcolor)
    text(1.05,-.8,'\leftarrowAvg. Cell Boundary','fontsize',16,'color',fcolor)
    
    set(gca,'fontsize',28,'XColor',fcolor,'YColor',fcolor,'YMinorTick','on','LineWidth',2)
    xt = 'Location on Trace (AU)';% input('enter the xaxis label','s');
    yt = 'Displacement (\mum)'; %input('enter the yaxis label','s');
    tt = 'Line-Profile Displacements';%input('enter the title','s');
    le = 'Shear'; %input('enter the legend','s');
    le2 = 'Normal';
    le3 = 'Border';
    xl = xlabel(xt);
    yl = ylabel(yt);
    %tl = title(tt);
    ytickformat('%.1f')
    set(xl, 'fontweight','bold','fontsize',28,'color',fcolor);
    set(yl,'fontweight','bold','fontsize',28,'color',fcolor);
    leg = legend([p1 p2],['\color{' fcolor '}' le],['\color{' fcolor '}' le2],'location','northwest');
    leg.FontSize = 20;
    legend boxoff
    %set(tl,'fontweight','bold','fontsize',title_font_size)
    axis([0 2 min(min(min(profBook2))) max(max(max(profBook2)))])
    
    
    title = strcat(prefix,['All_Data ' fcolor ' on ' bcolor]);
    savefile = [ListPath title];
    export_fig(ProfileOverlays,savefile,'-native');
    
    
    % Normalized Graph
    profBook3 = profBook2;
    for i = 1:size(profBook2,3)
        top(i) = abs(max(max(max(profBook2(2:3,:,i)))));
        profBook3(2:3,:,i) = profBook3(2:3,:,i)/top(i);
    end
    ProfileOverlaysNorm = figure;
    hold on
    for i = 1:size(profBook3,3)
        plot(profBook3(1,:,sortIdx(i))*aScale,profBook3(2,:,sortIdx(i)),'color',[map(i,1:3)])
        plot(profBook3(1,:,sortIdx(i))*aScale,profBook3(3,:,sortIdx(i)),'color',[map(i,1:3)],'linestyle','-.')
    end
    p1 = plot(profBook3(1,:,sortIdx(i))*aScale,mean(profBook3(2,:,sortIdx(:)),3,'omitnan'),'color',fcolor,'linewidth',3);
    p2 = plot(profBook3(1,:,sortIdx(i))*aScale,mean(profBook3(3,:,sortIdx(:)),3,'omitnan'),'color',fcolor,'linestyle','-.','linewidth',3);
    p3 = plot([cbf2 cbf2]*aScale,[min(min(min(profBook3))),1],'color',[.5 .5 .5],'linestyle','--','linewidth',1);
    p4=plot([.4 .4],[min(min(min(profBook2))),.6],'color',[1,.5,.5],'linestyle','--','linewidth',1);
    p5=plot([1.6 1.6],[min(min(min(profBook2))),max(max(max(profBook2)))],'color',[1,.5,.5],'linestyle','--','linewidth',1);
    set(gca,'Color',bcolor)
    text(1.05,-.25,'\leftarrowAvg. Cell Boundary','fontsize',16,'color',fcolor)
    
    set(gca,'fontsize',28,'XColor',fcolor,'YColor',fcolor,'YMinorTick','on','LineWidth',2)
    xt = 'Location on Trace (AU)';% input('enter the xaxis label','s');
    yt = 'Displacement (AU)'; %input('enter the yaxis label','s');
    tt = 'Line-Profile Displacements';%input('enter the title','s');
    le = 'Shear'; %input('enter the legend','s');
    le2 = 'Normal';
    le3 = 'Border';
    xl = xlabel(xt);
    yl = ylabel(yt);
    %tl = title(tt);
    set(xl, 'fontweight','bold','fontsize',28,'color',fcolor);
    set(yl,'fontweight','bold','fontsize',28,'color',fcolor);
    ytickformat('%.2f')
    leg = legend([p1 p2],['\color{' fcolor '}' le],['\color{' fcolor '}' le2],'location','northwest');
    leg.FontSize = 20;
    legend boxoff
    axis([0 2 min(min(min(profBook3))) max(max(max(profBook3)))])
    ytickformat('%.1f')
    %set(tl,'fontweight','bold','fontsize',title_font_size)
    
    
    title = strcat(prefix,['All_Data_Normalized ' fcolor ' on ' bcolor]);
    savefile = [ListPath title];
    export_fig(ProfileOverlaysNorm,savefile,'-native');
end