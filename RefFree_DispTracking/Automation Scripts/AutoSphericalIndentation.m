%%Automating Spherical Indentation Scripts
%% For Running Spherical Indentation Script
for i = 1:size(dirList,1)
    try
        SphericalIndent(dirList{i,1})
        disp(num2str(i))
    catch
        disp('Indentation Failed')
        i
    end
end

%%
for i = 1:size(dirList,1)
    cd(dirList{i})
    load('SphereIndent.mat')
    figure
    % Heatmap for XZ
    colorMapRad = colormap(jet(77000));
    colorMapRad(65537:end,:) = [];
    colorMapRad = flipud(colorMapRad);
    MaximumHeatMap = imagesc(radXX(:,1),radYY(:,1),radvq);
    radHeat = MaximumHeatMap.CData;%.*(imageBinary==0);
    
    radHeatNaN = (isnan(radHeat));
    radHeat(isnan(radHeat)) = 0;
    heatScale = (65536/2.4); %(max(max(imageHeat)))
    radHeat = uint16(round(radHeat * heatScale));
    radHeatColor = ind2rgb(radHeat,colorMapRad);
    
    
    scale = (max(max(radYY))-min(min(radYY)))/51;
    addTop = round(min(min(radYY)) / scale);
    addbottom = round((22 - max(max(radYY)))/ scale);
    total = zeros(addTop+addbottom+51,500,3);
    total(addTop:addTop+50,:,:) = radHeatColor(:,:,:);
    
    radHeatColor2 = imresize(total, [100 500]);
    
    radHeatColor3(:,:,:,i) =  radHeatColor2;
    close
    % figure
    % imshow(radHeatColor)
    % hold on
    % text(200,25,[num2str(max(max(radYY))) ' ' num2str(min(min(radYY)))])
    
    figure
    imshow(radHeatColor2)
    hold on
    text(200,25,[num2str(max(max(radYY))) ' ' num2str(min(min(radYY)))])
    
    
end
radHeatColor3(radHeatColor3(:,:,1:3,:)==0) = NaN;
rHC3Mean = mean(radHeatColor3(:,:,:,:),4,'omitnan');
mInd = figure;
imshow(rHC3Mean)
hold on
for i = 1:4
    plot([-1 525],[(500/22)*(i) (500/22)*(i)],'k')
    %text(200,25,[num2str(max(max(radYY))) ' ' num2str(min(min(radYY)))])
end
for i = 1:11
    plot([(2500/50)*(i) (2500/50)*(i)],[-1 105],'k')
    %text(200,25,[num2str(max(max(radYY))) ' ' num2str(min(min(radYY)))])
end
title = ['\Mean Indentation Profile'];
savefile = [ListPath title];
export_fig(mInd,savefile,'-native');
cd(ListPath)
