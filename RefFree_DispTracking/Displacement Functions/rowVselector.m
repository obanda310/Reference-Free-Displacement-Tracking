function [rowV] = rowVselector(image)

prompt = msgbox('Unable to automatically determine dot orientation! Please draw a line segment in the "horizontal" direction connecting at least two marker centers. When finished, double-click the line segment to save the line.')
uiwait(prompt)
rowVf = figure;
imshow(image.Proj,[]);
h = imline;
wait(h);
rowLine = getPosition(h);
rowV(1,[1 2]) = rowLine(1,1:2)-rowLine(2,1:2);
rowV(1,1) = rowV(1,1)*-1
close

end