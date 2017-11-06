%% Automated Script
clear all
close all
[ListName,ListPath] = uigetfile;
fullListName = strcat(ListPath,ListName);
load(strcat(ListPath,ListName),'dirList')
parts = strsplit(ListName, 'List');
prefix = parts{1};
 set(0,'defaultfigurecolor',[1 1 1])
%% For Running XY Script
for i = 1:size(dirList,1)
    dispXYfunc(dirList{i,1})
    disp(num2str(i))
end

%% For Running Z Script
for i = 1:size(dirList,1)
    cd(dirList{i,1})  
    %First check to see if the stack to analyze has been stored previously
    clear files
    files = dir('*.mat');
    if size(files,1)>=1
        for k = 1:length(files)
            current=files(k).name;
            if size(current,2)>12
            check(k)=strcmp(current(end-12:end),'StackName.mat');
            else
            check(k) = 0;
            end
        end
    end
    
    %If the Stack has not been assigned previously, manually assign it
    if size(find(check),2)==0
    [StackName,StackPath] = uigetfile('*.tif'); 
    save('StackName.mat','StackName')
    else
    load('StackName.mat','StackName')
    if StackName == 0
        [StackName,StackPath] = uigetfile('*.tif');
        save('StackName.mat','StackName')
    end
    end
    
end
%% Make sure every folder has a rowV file
[rowVName,rowVPath] = uigetfile;
load(strcat(rowVPath,rowVName),'rowV')
for i = 1:size(dirList,1)    
    cd(dirList{i,1})
    save('rowV.mat','rowV')
end

%% For Running Z Script
for i = 1:size(dirList,1)
    dispZfunc(dirList{i,1})
end

%% For Running Z Script Outputs Only (faster alternative)
for i = 1:size(dirList,1)
    dispZfuncOutputs(dirList{i,1})
end

%% For Running Profile Script

for i = 1:size(dirList,1)
try
[normXY,normZ,normAxis,cell_boundary] = profFunc(dirList{i,1});
profBook(1,1:size(normAxis,2),i) = normAxis;
profBook(2,1:size(normXY,2),i) = normXY;
profBook(3,1:size(normZ,2),i) = normZ;
cb(i) = cell_boundary;
save('Profile Data\ProfileOutputData.mat','normXY','normZ','normAxis')
catch
    caught(i) = 1;
    i
    dirList{i,1}
end
end
cbf = mean(cb);
save(strcat(ListPath,prefix,'Profiles.mat'),'profBook','cb','cbf')

%%

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

map = brewermap(size(profBook2,3),'*spectral');
ProfileOverlays = figure;
hold on
for i = 1:size(profBook2,3)
plot(profBook2(1,:,sortIdx(i))*aScale,profBook2(2,:,sortIdx(i)),'color',[map(i,1:3)])
plot(profBook2(1,:,sortIdx(i))*aScale,profBook2(3,:,sortIdx(i)),'color',[map(i,1:3)],'linestyle','-.')
end
p1= plot(profBook2(1,:,sortIdx(i))*aScale,mean(profBook2(2,:,sortIdx(:)),3,'omitnan'),'color',[0 0 0],'linewidth',3);
p2=plot(profBook2(1,:,sortIdx(i))*aScale,mean(profBook2(3,:,sortIdx(:)),3,'omitnan'),'color',[0 0 0],'linestyle','-.','linewidth',3);
p3=plot([cbf2 cbf2]*aScale,[min(min(min(profBook2))),max(max(max(profBook2)))],'color',[.3 .3 .3],'linestyle','--','linewidth',1);

text(1.05,-.8,'\leftarrowAvg. Cell Boundary','fontsize',16)

set(gca,'fontsize',28)
xt = 'Location on Trace(AU)';% input('enter the xaxis label','s');
yt = 'Displacement (\mum)'; %input('enter the yaxis label','s');
tt = 'Line-Profile Displacements';%input('enter the title','s');
le = 'Shear'; %input('enter the legend','s');
le2 = 'Normal';
le3 = 'Border';
xl = xlabel(xt);
yl = ylabel(yt); 
%tl = title(tt);

set(xl, 'fontweight','bold','fontsize',28); 
set(yl,'fontweight','bold','fontsize',28);
leg = legend([p1 p2],le,le2,'location','northwest');
leg.FontSize = 20;
%set(tl,'fontweight','bold','fontsize',title_font_size)
axis([0 2 min(min(min(profBook2))) max(max(max(profBook2)))])


title = strcat(prefix,'All_Data');
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
p1 = plot(profBook3(1,:,sortIdx(i))*aScale,mean(profBook3(2,:,sortIdx(:)),3,'omitnan'),'color',[0 0 0],'linewidth',3);
p2 = plot(profBook3(1,:,sortIdx(i))*aScale,mean(profBook3(3,:,sortIdx(:)),3,'omitnan'),'color',[0 0 0],'linestyle','-.','linewidth',3);
p3 = plot([cbf2 cbf2]*aScale,[min(min(min(profBook3))),1],'color',[.3 .3 .3],'linestyle','--','linewidth',1);

text(1.05,-.3,'\leftarrowAvg. Cell Boundary','fontsize',16)

set(gca,'fontsize',28)
xt = 'Location on Trace(AU)';% input('enter the xaxis label','s');
yt = 'Displacement (AU)'; %input('enter the yaxis label','s');
tt = 'Line-Profile Displacements';%input('enter the title','s');
le = 'Shear'; %input('enter the legend','s');
le2 = 'Normal';
le3 = 'Border';
xl = xlabel(xt);
yl = ylabel(yt); 
%tl = title(tt);
set(xl, 'fontweight','bold','fontsize',28); 
set(yl,'fontweight','bold','fontsize',28);
leg = legend([p1 p2],le,le2,'location','northwest');
leg.FontSize = 20;
axis([0 2 min(min(min(profBook3))) max(max(max(profBook3)))])

%set(tl,'fontweight','bold','fontsize',title_font_size)


title = strcat(prefix,'All_Data_Normalized');
savefile = [ListPath title];
export_fig(ProfileOverlaysNorm,savefile,'-native');

