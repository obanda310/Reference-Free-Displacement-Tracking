[file,path] = uigetfile('*.tif');
filepath = [path,file];
% filepath = 'C:\Users\Janty\Desktop\OmarStacks\20150703  Z stack 633+transmitted dot array 20%laser 18um depth_z003_c001.tif';
img = imread(filepath);
figure
imshow(img,[]) % Select a "c001" file
% Eventually, should loop through each image in stack, i.e. z001, zoo2,
% etc.
%%
% Make binary image and fill in gaps
% level = graythresh(img);
level = graythresh(img);
bw = im2bw(img,level);
bw = imfill(bw,'holes');
figure
imshow(bw,[])