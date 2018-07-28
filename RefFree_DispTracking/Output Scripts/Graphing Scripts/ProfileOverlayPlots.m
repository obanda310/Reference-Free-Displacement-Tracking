clear all
close all

%addtopath('Export_Fig')

%CHANGE FONT SIZES HERE
AxisFontSize = 28;
AxisTitleFontSize = 28;
LegendFontSize = 20;




prefix = '/5%wt Gels ';
ListPath = cd;
load('5%wt GelsProfiles.mat')


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


set(gca,'fontsize',AxisFontSize,'LineWidth',2)
xt = 'Location on Trace(AU)';% input('enter the xaxis label','s');
yt = 'Displacement (\mum)'; %input('enter the yaxis label','s');
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

set(gca,'fontsize',AxisFontSize)
xt = 'Location on Trace(AU)';% input('enter the xaxis label','s');
yt = 'Displacement (AU)'; %input('enter the yaxis label','s');
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
axis([0 2 min(min(min(profBook3))) max(max(max(profBook3)))])

%set(tl,'fontweight','bold','fontsize',title_font_size)


title = strcat(prefix,'All_Data_Normalized');
savefile = [ListPath title];
export_fig(ProfileOverlaysNorm,savefile,'-native');