function [heatmap] = Auxheatmap(xqMax,yqMax,vqIni,LUT,filename,savePath,maxD,CO,filtMask)
if nargin == 6
maxD = max(abs(vq(:)));
end
xq = 1:1:xqMax;
yq = 1:1:yqMax;
vq = imresize(vqIni,[xqMax yqMax]);
vqfilt = imresize(filtMask,[xqMax yqMax]);
vq(abs(vq)<CO & filtMask) = 0;
%Check to see if data includes direction or is simply magnitudes
sidedCheck = min(vqIni(:))<0;

%Build heat map accordingly
if sidedCheck == 1
    %doublesided
    scaleD = 32768/maxD;
    map = brewermap(65535,LUT);
    hm = imagesc(xq,yq,vq);
    hm2 = hm.CData;
    hm2(hm2>0) = 32768+(abs(hm2(hm2>0))*scaleD);
    hm2(hm2<0) = 32768 - (abs(hm2(hm2<0))*scaleD);
    hm2(hm2==0) = 32768;
    hm2(isnan(hm2)) = 32768;
    hm2 = uint16(hm2);
else
    %singlesided
    scaleD = 65535/maxD;
    map = brewermap(65535,LUT);
%     if strcmp(LUT,'*spectral')==1
%       map(1,:) = [0,0,0];  
%     end
    hm = imagesc(xq,yq,vq);
    hm2 = hm.CData;
    hm2 =(hm2*scaleD);
    hm2(isnan(hm2)) = 0;
    hm2 = uint16(hm2);
end
imageHeatColor = single(ind2rgb(hm2,map));
heatmap(:,:,:) = imageHeatColor;
close
maxHeatMap = figure;
hold on
imshow(imageHeatColor);

%% Export Image
savefile = [savePath filename '.tif'];
export_fig(maxHeatMap,savefile,'-native');

%% Colorbar
mkdir('HeatMaps\Traction','ColorBar')


colorBar1 = single(zeros(500,25));
range = uint16(round(linspace(65536,1,500)'));
for i = 1:25
    colorBar1(1:500,i) = range;
end
colorBar2 = ind2rgb(colorBar1,map);
for i = 1:10
    colorBar2((i*50)-3:(i*50),13:25,:) = 0;
end

colorBar2(1:3,13:25,:) = 0;
%Save Color Bar Image

colorBarSave = figure;
hold on
imshow(colorBar2);
savefile = [savePath '\ColorBar\' filename 'ColorBar.tif'];
export_fig(colorBarSave,savefile,'-native');

close

filePath=cd;
cd HeatMaps\Traction\ColorBar
colorBarTxt = fopen([filename 'ColorBarValues.txt'],'wt');
p1Format = 'Tick number %1.0f is %.2f \n';
if sidedCheck == 1
    for i = 1:11
        colorBarValues(i,1) =  round(maxD - ((maxD/5)*(11-i)),2); % round(max(vq(:)) - ((max(vq(:))/10)*(11-i)),2)
        fprintf(colorBarTxt,p1Format,i,colorBarValues(i,1));
    end
else
    for i = 1:11
        colorBarValues(i,1) =  round(maxD - ((maxD/10)*(11-i)),2); % round(max(vq(:)) - ((max(vq(:))/10)*(11-i)),2)
        fprintf(colorBarTxt,p1Format,i,colorBarValues(i,1));
    end
    
end
fclose(colorBarTxt);
cd(filePath)
