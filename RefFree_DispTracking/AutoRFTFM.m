%% Automated Script
%% Grab Directory List
clear all
close all
[ListName,ListPath] = uigetfile;
fullListName = strcat(ListPath,ListName);
load(strcat(ListPath,ListName),'dirList')
parts = strsplit(ListName, 'List');
prefix = parts{1};
set(0,'defaultfigurecolor',[1 1 1])

%% Clean Old Files in Directories
clear failed
for i = 1:size(dirList,1)
    try
        MoveMatFiles(dirList{i,1})
        disp(num2str(i))
    catch        
        failed{i} = dirList{i,1};
    end
end

for i = 1:size(dirList,1)
    try
        clearDirectory(dirList{i,1})
        disp(num2str(i))
    catch        
        failed{i} = dirList{i,1};
    end
end


%% Reference-Free Deformation Estimation Script
clear failed
for i = 1:size(dirList,1)
     try
        RFTFM(dirList{i,1})
        disp(num2str(i))
     catch
         disp('RFTFM Failed')
         failed(i) = i;
    end
end

%% Update Deformation Outputs (for reformatting data)
clear failed
for i = 1:size(dirList,1)
    try
        UpdateOutputs(dirList{i,1})
        disp(num2str(i))
    catch
        disp('RFTFM Failed')
        failed(i) = i;
    end
end
%% Tractions

clear failed
progressbar([prefix ' Calculating Tractions'])
for i = 1:size(dirList,1)
    progressbar(i/size(dirList,1))
    try
        [output] = interpFinal3D_2(dirList{i,1});
        disp(num2str(i))
    catch
        disp('RFTFM Failed')
        failed{i} = dirList{i,1};
    end
end
%% Cropping Images

%% Splitting Channels 1
clear failedSplit
for i = 1:size(dirList,1)
    try
        cd(dirList{i,1})
        %Set directory to save output split channels
        OutputDirectory = dirList{i,1};
        
        %set any changes in directory to find raw files here
        cd('Raw Images')
        FileDirectory = cd;
        
        %Name each channel here
        ch1 = 'Raw Transmitted Offset';
        ch2 = 'Raw Talin-GFP Offset';
        ch3 = 'Raw pFAK AF405 Offset';
        identfierSuffix = 'AF405-1.tif';        
        SplitChannels(FileDirectory,OutputDirectory,identfierSuffix,ch1,ch2,ch3)
    catch
        failedSplit(i) = i;
    end
end
%% Splitting Channels 2
clear failedSplit
for i = 1:size(dirList,1)
    try
        cd(dirList{i,1})
        %Set directory to save output split channels
        OutputDirectory = dirList{i,1};
        
        %set any changes in directory to find raw files here
        cd('Raw Images')
        FileDirectory = cd;
        
        %Name each channel here
        ch1 = 'Empty Channel';
        ch2 = 'Raw Transmitted True';
        ch3 = 'Raw Talin-GFP True';
        
        identfierSuffix = 'Transmitted.czi #2-1.tif';
        
        SplitChannels(FileDirectory,OutputDirectory,identfierSuffix,ch1,ch2,ch3)
    catch
        failedSplit(i) = i;
    end
end
%% Slice Picker
clear failedSlice
for i = 1:size(dirList,1)
    try
        keysuffix = 'FAK 576.tif';
        OutputDirectory = dirList{i,1};
        FileDirectory = dirList{i,1};
        
        SlicePicker(FileDirectory,OutputDirectory,keysuffix)
    catch
        failedSlice(i) = i;
    end
end
%% Determine transform needed to register images
clear failedTransform
close all
for i = 1:size(dirList,1)   
    try
        movekey = 'Raw Transmitted Offset_Slice';
        fixedkey = 'Raw Transmitted True_Slice';
        OutputDirectory = dirList{i,1};
        FileDirectory = dirList{i,1};
        OutputName = 'Raw Transmitted Offset_Registered.tif';
        [tform, moved] = registerImages(OutputDirectory,OutputName,FileDirectory,movekey,fixedkey);
    catch
        failedTransform(i) = i;
    end    
end

%% Analyzing Adhesions Pt 1
clear failedAna
%pFAK AF405 Offset_Slice
%FAK 576_Slice_
%Raw Talin-GFP Offset_Slice_

progressbar('Analyzing Adhesions')
for i = 1:size(dirList,1)
    progressbar(i/size(dirList,1))
    try
        cd(dirList{i,1});
        files = dir('*Raw Talin-GFP Offset_Slice_*');
        img = getImages(files(1).name);
        ProcessAdhesions(img,1,'Talin');
    catch
        failedAna{i,1} = dirList{i,1};
    end
end
%%
clear failedAna

progressbar('Analyzing Adhesions')
for i = 1:size(dirList,1)
    progressbar(i/size(dirList,1))
    try
        cd(dirList{i,1});
        files = dir('*pFAK AF405 Offset_Slice*');
        img = getImages(files(1).name);
        ProcessAdhesions(img,1,'FAK');
    catch
       failedAna{i,1} = dirList{i,1};
    end
end
%%
clear failedAna

progressbar('Analyzing Adhesions')
for i = 1:size(dirList,1)
    progressbar(i/size(dirList,1))
    try
        cd(dirList{i,1});
        files = dir('*FAK 576_Slice_*');
        img = getImages(files(1).name);
        ProcessAdhesions(img,0,'FAK');
    catch
        failedAna{i,1} = dirList{i,1};
    end
end

%% Analyzing Adhesions Pt 2
clear failedAna
for i = 1:size(dirList,1)
    try
        
        image1suffix = 'Raw pFAK AF405 Offset_Slice';
        image2suffix = 'Raw Talin-GFP Offset_Slice';        
        FileDirectory = dirList{i,1};
        [minmax] = grabMinMaxValues(FileDirectory,image1suffix,image2suffix);
        mmAll(i,1) = minmax(1,1) ;
        mmAll(i,2) = minmax(1,2);
    catch
        failedAna(i) = i;
    end
end
globalminmax(1,1) = min(mmAll(:,1));
mmAll(mmAll(:,2)==65535,:) = []; %remove entries that have 2^16 because they contain overexposed objects that probably are not adhesions
globalminmax(1,2) = mean(mmAll(:,2)) + 2*std(mmAll(:,2));
cd(ListPath)
save('GlobalMinMaxes.mat','globalminmax');
%%
clear failedAna
close all
for i = 1:size(dirList,1)   
%     try
        image1suffix = 'FAK';
        image2suffix = 'Talin';
        OutputDirectory = dirList{i,1};
        FileDirectory = dirList{i,1};
        [RatioMap,abBW] = RatiometricIFCompare(FileDirectory,OutputDirectory,image1suffix,image2suffix,[1 1]);        
%     catch
%         failedAna(i) = i;
%     end    
end
%%
clear failed
for i = 1:size(dirList,1)
%     try
        [IntVal,DefVal] = CompareIntensityDeformation(dirList{i,1});
        AllVals(1:size(IntVal,1),1,i) = IntVal;
        AllVals(1:size(IntVal,1),2,i) = DefVal;
        disp(num2str(i))
%     catch
%         disp('RFTFM Failed')
%         failed(i) = i;
%     end
end
%%
AllVals2 = [0,0];
for i = 1:size(AllVals,3)
    AllVals2 = cat(1,AllVals2,AllVals(:,:,i));
end
AllVals2(AllVals2(:,1)==0,:) = [];
AllVals2(AllVals2(:,2)<0.5,:) = [];
figure
scatter(AllVals2(:,1),AllVals2(:,2))

%% 

clear failed
largestSizeX = 0;
largestSizeY = 0;
for i = 1:size(dirList,1)
    try
        thisMap = CreateOverlay(dirList{i,1});
        thisSizeX = size(thisMap,1);
        thisSizeY = size(thisMap,2);
        if thisSizeX > largestSizeX
            largestSizeX = thisSizeX;
        end
        if thisSizeY > largestSizeY
            largestSizeY = thisSizeY;
        end
    catch       
        failed(i) = i;
    end
end
%%
for i = 1:size(dirList,1)
    try
        thisMap = CreateOverlay(dirList{i,1});
%         xStart = ; 
%         yStart = ;
    catch       
        failed(i) = i;
    end
end

%%
for i = 1:size(dirList,1)
    try
        getCellTrace(dirList{i,1})
    catch
        disp(['failed ' num2str(i)])
    end
end



%% Delete rowV.mat
for i = 1:size(dirList,1)    
   
   try
       cd(dirList{i,1})
       delete('Matlab Data Files\rowV.mat')
       disp('Success')
   catch
       disp('Failed')
   end
end

%% Grab Frame sizes
for i = 1:size(dirList,1)
    
    try
        cd(dirList{i,1})
        load('Matlab Data Files\3ddata.mat','image')
        imgSizes(i,1) = size(image.Area,1);
        imgSizes(i,2) = size(image.Area,2);
    catch
        disp(dirList{i,1})
    end
end
max(imgSizes)
%% Print and concatenate all Quiver-Adhesion Images

prefix = 'FAK'; %'Talin'; % 
for i = 1:size(dirList,1)
    try
        QuiverAdhesions(dirList{i,1},prefix)
    catch        
    end
end
%%
tic
prefix = 'Talin';
QAs = [];
Trans = [];
for i = 1:size(dirList,1)    
    try
        toc
        cd(dirList{i,1})
        text = dirList{i,1};
        temp = imread([prefix ' QuiverAdhesions.tif']);
        temp = insertText(temp,[10,10],text,'FontSize',20);
        QAs=cat(4,QAs,temp);
    catch    
        text = dirList{i,1};
        temp = uint8(zeros(1001,1001,3));
        temp = insertText(temp,[10,10],text,'FontSize',20);                
        QAs=cat(4,QAs,temp);
        
    end
end

for i = 1:size(dirList,1)    
    try
        toc
        cd(dirList{i,1})
        text = dirList{i,1}; 
        temp2 = imread('Extended Border Transmitted.tif');
        temp2 = insertText(temp2,[10,10],text,'FontSize',20);
        temp3 = temp2(:,:,1);
        Trans = cat(3,Trans,temp3);
    catch    
        text = dirList{i,1};               
        temp2 = uint8(zeros(1000,1000));
        temp2 = insertText(temp2,[10,10],text,'FontSize',20);
        temp3 = temp2(:,:,1);
        Trans = cat(3,Trans,temp3);
    end
end

ShowStackC(QAs)
cd(ListPath)
save([prefix 'QuiverAdhesions.mat'],'QAs')
for i = 1:size(QAs,4)
    thisImg(:,:,:) = QAs(:,:,:,i);
    imwrite(thisImg,[prefix 'QuiverAdhesions.tif'],'WriteMode','append');
    thisTImg(:,:) = Trans(:,:,i);
    imwrite(thisTImg,[prefix 'ExtTrans.tif'],'WriteMode','append');
end
%% Create all individual reports
TI = [];
TAA = [];
TA = [];
TDxy = [];
TDz = [];
TD3D = [];
TF = [];
TE3D = [];
TET = [];
TAAT = [];
TIT = [];

cd(ListPath)
mkdir('Graphs and Images')
cd('Graphs and Images')
outputDirectory = cd;  
for i = 1:size(dirList,1)
     try
        [Total] = CreateAdhesionReport(dirList{i,1},ListPath,outputDirectory,i);
        TI = cat(1,TI,Total.Int);
        TAA = cat(1,TAA,Total.AdhArea);
        TA = cat(1,TA,Total.SpreadArea);
        TDxy = cat(1,TDxy,Total.DeformationXY);
        TDz = cat(1,TDz,Total.DeformationZ);
        TD3D = cat(1,TD3D,Total.Deformation3D);
        TF = cat(1,TF,Total.Force);
        TE3D = cat(1,TE3D,Total.Energy3D);
        TET = cat(1,TET,Total.EnergyTop);
        if Total.td ==1
            td = 1;
            TAAT = cat(1,TAAT,Total.AdhAreaTalin);
            TIT = cat(1,TIT,Total.IntTalin);
        else
            td = 0;
        end
     catch
     end
end

cd(ListPath)
if td == 1
save([prefix 'StatTotals.mat'],'TI','TA','TAA','TDxy','TDz','TD3D','TF','TE3D','TET','TIT','TAAT','td')
else
save([prefix 'StatTotals.mat'],'TI','TA','TAA','TDxy','TDz','TD3D','TF','TE3D','TET','td')
end
%% Create Summary Reports


