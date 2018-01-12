function ColorScheme(fcolor,bcolor,label,le,AxisFontSize,LegendFontSize,box,AzEle)
if AzEle(1,1) ~=360
    view([AzEle(1,1) AzEle(1,2)])
end
set(gcf,'color',bcolor)
set(gca,'Color',bcolor)
set(gca,'fontsize',AxisFontSize,'XColor',fcolor,'YColor',fcolor,'ZColor',fcolor)

for j = 1:size(label,2)
    set(label{j}, 'fontweight','bold','fontsize',28,'color',fcolor);
end
if strcmp(le{1},'0')~=1
    for j = 1:size(le,2)
        leg{j} = ['\color{' fcolor '}' le{1,j}];
    end
    
    legend(leg,'location','best','FontSize',14);
    %leg.FontSize = LegendFontSize;
    if box ~=1
        legend boxoff
    end
    
end
%set(tl,'fontweight','bold','fontsize',title_font_size)

end