function p = pub_fig(x,y,line_width,line_style,line_color,marker,marker_size, marker_edge_color,marker_face_color,title)

ProfLinePlotClean = figure;
hold on
axis_label_font_size = 20;
title_font_size = 24;

for i = 1:size(x,2)
p = line(x(:,i),y(:,i));
set(p,'linewidth',line_width,'linestyle',line_style,'color',line_color{i,1});
%set(p,'marker',marker,'markeredgecolor',marker_edge_color{i,1},'markerfacecolor',marker_face_color{i,1},'markersize',marker_size);
end
set(gca,'fontsize',16)
xt = 'Location on Trace(microns)';% input('enter the xaxis label','s');
yt = 'Displacement (microns)'; %input('enter the yaxis label','s');
tt = 'Line-Profile Displacements';%input('enter the title','s');
le = 'Shear Trace'; %input('enter the legend','s');
le2 = 'Normal Trace';
xl = xlabel(xt);
yl = ylabel(yt); 
%tl = title(tt);
leg = legend(le,le2);
set(xl, 'fontweight','bold','fontsize',axis_label_font_size); 
set(yl,'fontweight','bold','fontsize',axis_label_font_size);
%set(tl,'fontweight','bold','fontsize',title_font_size)

filePath = cd;
savefile = [filePath strcat('\Profile Data\',title)];
export_fig(ProfLinePlotClean,savefile,'-native');
