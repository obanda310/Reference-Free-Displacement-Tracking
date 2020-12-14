function A6 = removeLarge2(A,scale)
%% Minimum cutoff
A2 = A;
A2(A2<prctile(A(:),10))=0;
%% Threshold Minimum Cutoff
A3 = A2>0;
%% Erode Small Features
A4 = A3;
SE = strel('disk',3);

for i = 1:7
A4 = imerode(A4,SE);
end

for i = 1:7
A4 = imdilate(A4,SE);
end
%%
SE = strel('disk',3,0);
A4p = A3;
for i = 1:size(A,3)
%     A4p(4:7,:,i) = 1;
%     A4p(end-6:end-3,:,i) = 1;
%     
%     A4p(1:3,:,i) = 0;
%     A4p(end-2:end,:,i) = 0;

    A4p(:,:,i) = imerode(imfill(imdilate(bwmorph(edge(A4p(:,:,i),'Canny',[.1 .2],15),'close'),SE),'holes'),SE);
    
    A4p(7:end-7,1:3,i) = 1;
    A4p(7:end-7,end-2:end,i) = 1;
    %A4p(1:3,7:end-7,i) = 1;
    %A4p(end-2:end,7:end-7,i) = 1;
    
     A4p(:,:,i) = imerode(imfill(imdilate(A4p(:,:,i),SE),'holes'),SE);
     A4p(:,:,i) = imdilate(imerode(A4p(:,:,i),SE),SE);
end

%%
A4t = A4p|A4;


%% Use A4 mask to filter A

A5 = A.* (A4t==0);
%% Resize image to standard 2048x2048 size that I have been using for everything

A6 = imresize(A5,[size(A5,1)*scale size(A5,2)*scale], 'bicubic');
    
