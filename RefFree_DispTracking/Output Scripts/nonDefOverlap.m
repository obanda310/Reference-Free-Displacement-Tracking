function percentcovered = nonDefOverlap(directory)
if nargin == 1
    cd(directory);
end
%%
clear overlap
load('3Ddata.mat')
vq3 = double(vq3);
vq3(isnan(vq3)) = 0;
maskvq3 = double((image.ADil)==0);
for i = 1:size(vq3,3)   
overlap(:,:,i) = abs(vq3(:,:,i)).*maskvq3;
end
sumpre = sum(sum(sum(abs(vq3))));
sumpost = sum(sum(sum(overlap)));
percentcovered = sumpost/sumpre
ShowStack(overlap)
