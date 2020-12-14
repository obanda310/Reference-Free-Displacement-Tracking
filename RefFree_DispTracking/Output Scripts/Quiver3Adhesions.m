function Quiver3Adhesions(directory)
if nargin == 1
    cd(directory)
end
%%
close all
load('Matlab Data Files\3Ddata.mat')
load('Matlab Data Files\FAK Adhesions.mat')

Xlimits = [r.s(1,2),r.s(1,2);0,0];
Ylimits = [0,r.s(2,2);0,r.s(2,2)];
Zlimits = Surface2(Ylimits,Xlimits);
%%
figure
quiver3(r.rX,r.rY,r.rZ,r.dX,r.dY,r.dZ)
hold on

h = surface(Ylimits,Xlimits,Zlimits,imgO, 'facecolor', 'texturemap','facealpha','texturemap', 'edgecolor', 'none');
h.AlphaData = imgO/prctile(imgO(:),90);
%%
figure
q = quiver3(r.rX,r.rY,r.rZ,r.dX,r.dY,r.dZ);
hold on
image.MaskStack(image.MaskStack<300) = 0;
for i = 1:size(image.MaskStack,3)
    Xlimits = [r.s(1,2),r.s(1,2);0,0];
    Ylimits = [0,r.s(2,2);0,r.s(2,2)];
    Zlimits = (0.4*i)*ones(2,2);
    h = surface(Ylimits,Xlimits,Zlimits,image.MaskStack(:,:,i), 'facecolor', 'texturemap','facealpha','texturemap', 'edgecolor', 'none');
    h.AlphaData = 0.1*(image.MaskStack(:,:,i)>0);
   
end
%%
img2 = img;
img2(img2<1000) = 0;
Xlimits = [r.s(1,2),r.s(1,2);0,0];
Ylimits = [0,r.s(2,2);0,r.s(2,2)];
Zlimits = Surface2(Ylimits,Xlimits);
h2 = surface(Ylimits,Xlimits,Zlimits,img2, 'facecolor', 'texturemap','facealpha','texturemap', 'edgecolor', 'none');
h2.AlphaData = 0.7*(img2>0);