function polygons = bwperimtrace_frontend(im)

im = im>0;
im = bwlabel(im);

xsize = 700;
ysize = 700;

xv = (1:xsize)-(xsize/2);
yv = (1:ysize)-(ysize/2);

bounds = cell(0);

indx = unique(im(:))';
indx(indx==0)=[];

for t = indx
    t
    bounds{end+1} = bwperimtrace(im == t,[xv(1) xv(end)],[yv(1) yv(end)]);
end

%%
clear polygons
for i = 1:100
   polygons{i,1} = bounds{1,i}{1,1}(:,:); 
end 
