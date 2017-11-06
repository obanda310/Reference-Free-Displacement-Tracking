%% Making the spherical indentation data for FE fitting
h = fspecial('average', [20 20]);
vqtest(isnan(vqtest)) = 0;
for i = 1:size(vqtest,3)
vqmask(:,:,i) = bwareaopen(vqtest(:,:,i)<0,50000);
end
vqtestfilt1 = vqtest.*double(vqmask);
for i = 1:size(vqtest,3)
    vqtestfilt2(:,:,i) = filter2(h, vqtestfilt1(:,:,i));
end
vqtestfilt3 = vqtestfilt2*-1;
ShowStack(vqtestfilt3)
%%
vqtestproj = squeeze(mean(permute(vqtestfilt3,[3,1,2])));
imshow(vqtestproj,[])

%inputs:
%originalmatrix: the original matrix
%binarisedmatrix = originalmatrix > threshold; %the thresholded matrix, a logical array
[rows, cols] = ndgrid(1:size(vqtestproj, 1), 1:size(vqtestproj, 2));
rc = sum(rows(vqtestproj>0) .* vqtestproj(vqtestproj>0)) / sum(vqtestproj(vqtestproj>0));
cc = sum(cols(vqtestproj>0) .* vqtestproj(vqtestproj>0)) / sum(vqtestproj(vqtestproj>0));

hold on
scatter(cc,rc,50)

%%
for j = 1:size(vqtestfilt3,3)
for i = 1:500   
vqmask2 = double(not(vqtestfilt3(:,:,j)>=0));
vqmask2 = rgb2gray(insertShape(vqmask2,'circle',[cc rc i],'LineWidth',1))>0;
vqmask2area = sum(sum(single(vqmask2>0)));
radprofile(i,j) = sum(sum(vqmask2.*vqtestfilt3(:,:,j)))/vqmask2area;
i
    end
end
figure
hold on
%%
for i = 1:size(vqtestfilt3,3)
plot(radprofile(:,i))
end

%%
planesLocs3 = [6.53,11.27,15.90];
for i = 1:size(vqtestfilt3,3)
    radprofile2(:,i) = radprofile(:,i) + planesLocs3(i);
end
radprofile3 = radprofile(:,1:3);
radprofileXs = [[1:1:500]' * xyScale] *ones(1,3);
radprofileZs = ones(500,1) * planesLocs3  ;

[radXX radYY] = meshgrid(([1:1:500]' * xyScale), min(planesLocs3):(max(planesLocs3)-min(planesLocs3))/50:max(planesLocs3));
radvq = griddata(radprofileXs,radprofileZs,radprofile3,radXX,radYY);
%% Heatmap for XZ
colorMapRad = brewermap(65536,'spectral');
MaximumHeatMap = imagesc(radXX(:,1),radYY(:,1),radvq);
radHeat = MaximumHeatMap.CData;%.*(imageBinary==0);

radHeatNaN = (isnan(radHeat));
radHeat(isnan(radHeat)) = 0;
heatScale = (65536/max(radprofile3(:))); %(max(max(imageHeat)))
radHeat = uint16(round(radHeat * heatScale));
radHeatColor = ind2rgb(radHeat,colorMapRad);
%%
close all
imshow(radHeatColor)
hold on 
plot([5 10 15 20 25 30 35 40 45 50 55 60 65 70]/xyScale,[1 50 1 50 1 50 1 50 1 50 1 50 1 50] )