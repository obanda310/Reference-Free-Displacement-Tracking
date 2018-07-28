function savestack(Stack)
Stack=uint16(Stack);
[StackName,StackPath] = uiputfile('*.tif');
for i = 1:size(Stack,3)
imwrite(Stack(:,:,i),StackName,'WriteMode','append');                   
end