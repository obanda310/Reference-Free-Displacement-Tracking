stack = a.images;
img = stack(:,:,20);
imshow(img,[])
r = imrect(gca);
rPos = round(getPosition(r));

noImgs = size(stack,3);
newStack = zeros(rPos(4)+1,rPos(3)+1,noImgs);
for i = 1:noImgs
    newImg = imcrop(stack(:,:,i),rPos);
    newStack(:,:,i) = newImg;
end
% f = msgbox('Double-Click to crop after making rectangular selection!');
% waitfor(f);
% handles = guidata(gcf);
% OutputStack = getappdata(hFig,'MyMatrix');
% CurrentFrame = round((get(handles.SliderFrame,'Value')));
% set(handles.Edit1,'String',num2str(CurrentFrame));    
% 
% [~, bounds] = imcrop(OutputStack(:,:,CurrentFrame));
% bounds = round(bounds);
% CroppedStack = OutputStack(bounds(1,2):bounds(1,2)+bounds(1,4),bounds(1,1):bounds(1,1)+bounds(1,3),:);
% cropImages = Original(bounds(1,2):bounds(1,2)+bounds(1,4),bounds(1,1):bounds(1,1)+bounds(1,3),:);
% setappdata(hFig,'MyMatrix',CroppedStack);
% setappdata(hFig,'MyOrigMatrix',cropImages);
% OutputStack = getappdata(hFig,'MyMatrix');
% imshow(OutputStack(:,:,CurrentFrame),[],'Parent',handles.axes1);
% guidata(hFig,handles);