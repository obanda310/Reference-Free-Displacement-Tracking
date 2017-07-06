close all
for i = 1:20
imagesBpass(:,:,i) = bpass(image2,1,i*5)>0;
end

ShowStack(imagesBpass)