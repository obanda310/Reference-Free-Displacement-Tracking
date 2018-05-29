function [Ys, Xs, Ys2, Xs2, rad]=RadialProfiles(directory,type)
if nargin == 1 || nargin == 2
    cd(directory);
    if strcmp(type,'pill ') == 1
        type = 'pill';
    elseif strcmp(type,'circle ') == 1
        type = 'circle';
    end
elseif nargin == 0
    type = input('pill or circle?','s');
end
%%

close all
filePath = cd;
load('Shear Mat Files\ShearCutout.mat')
load('Shear Mat Files\ShearImages.mat')

%%
if strcmp(type,'pill') == 1
    %% Create a trace of a pill shape end protuding from left side
    D = bwdist(image.Area); %Dist Transform
    E(:,:) = D(:,:)>prctile(D(:),99); %Find 'center line'
    
    %Center of Semi-Circle End
    [~,ERx] = find(E,1,'last');
    E2y = find(E(:,ERx)>0);
    ERy = mean(E2y);
    pts(1,1) = ERx;
    pts(1,2) = ERy;
    
    %Horizontal Edge for approximating 0 degrees
    ELx = 1;
    ELy = mean(find(E(:,ELx)>0));
    
    %Length of Trace micronse
    lTrace = 30; 
    
    
    
    %% For the circle case
elseif strcmp(type,'circle') == 1
    %%
    clear Di Db D
    %Need to define: pts, pts2, ERy, ERx, ELy, ELx
    Di = bwlabel(image.Area==0);
    for i = 1:max(Di(:))
        sizeCounts = size(find(Di==i),1);
    end
    D = Di(:,:)==find(sizeCounts == max(sizeCounts),1,'first');
    Db = bwboundaries(D);
    for i = 1:size(Db{1},1)
        clear distances
        distances = (Db{1}(:,1) - Db{1}(i,1)).^2 + (Db{1}(:,2) - Db{1}(i,2)).^2;
        maxDistance(i,1) = max(distances);
        maxDistance(i,2) = i;
        maxDistance(i,3) = find(distances == max(distances),1,'first');
    end
    
    idx = find(maxDistance(:,1) == max(maxDistance(:,1)),1,'first');
    pts(1,1) = mean( [Db{1}(maxDistance(idx,2),2) Db{1}(maxDistance(idx,3),2)]);
    pts(1,2) = mean( [Db{1}(maxDistance(idx,2),1) Db{1}(maxDistance(idx,3),1)]);
    
    %Center of Circle
    ERx = pts(1,1);
    ERy = pts(1,2);
    
    %Horizontal Edge for approximating 0 degrees
    ELx = size(image.Area,2);
    ELy = ERy;
    
    %Length of Trace microns
    lTrace = 50;
    
    %Determine the radius of the circle region (assume masDistance is
    %diameter)
    Crad = sqrt(maxDistance(idx,1))/2;
    
    
end
lTP = lTrace/.1625; %Lenth of trace in pixels

%% Create Thick Pie Shaped Profiles Through Max Shear
%Code should get line profiles for several lines parallel to the first line
%and average them for "potentially" cleaner data
clear dispProfXYThck

profWidth = 20; %In Degrees

%Coordinates of Max Shear
[SMy,SMx] = find(vq3 == max(max(vq3)));
pts(2,1) = SMx;
pts(2,2) = SMy;

%Coordinates of Last point on Trace through Max Shear. Found by scaling
%vector from center to max shear so that Length(Profile) = lTP.
dx = pts(2,1)-pts(1,1);
dy = pts(2,2)-pts(1,2);
lxy = sqrt(dx^2+dy^2);
lxyScale = lTP/lxy;
pts(2,1) = pts(1,1) + lxyScale*dx;
pts(2,2) = pts(1,2) + lxyScale*dy;

% Plot Relevant Points and Profiles on Binary Image
figure
imshow(image.Area)
hold on
scatter(pts(1,1),pts(1,2))
scatter(ELx,ELy)
scatter(SMx,SMy)

for i = 1:profWidth
    %Iteritively increment by 1 degree and take profile
    ptsTemp = pts;
    ptsTemp(2,1) = pts(1,1)+(pts(2,1)-pts(1,1))*cos((pi/180)*(i-profWidth/2))+(pts(2,2)-pts(1,2))*sin((pi/180)*(i-profWidth/2));
    ptsTemp(2,2) = pts(1,2)-(pts(2,1)-pts(1,1))*sin((pi/180)*(i-profWidth/2))+(pts(2,2)-pts(1,2))*cos((pi/180)*(i-profWidth/2));
    clear profTemp
    profTemp = improfile(vq3,ptsTemp(:,1),ptsTemp(:,2),'bilinear');
    profTemp(profTemp==0,:) = []; %Removes Zeros to start profile right at shape edge
    dispProfXYThck(1:size(profTemp,1),i)= profTemp(:,1);
    plot(ptsTemp(1:2,1),ptsTemp(1:2,2))
    
end
dispProfXYThckFinal = mean(dispProfXYThck,2);


Ys = [dispProfXYThckFinal];

%% Plot Max Profile
Xs = 0:.1625:(.1625*(size(Ys,1)-1));
figure
plot(Xs,Ys)

%% Do it over again. This time sample at increments starting from -90degrees

%Determine Length of Profile Trace
Edx = ERx - ELx;
Edy = ERy - ELy;
Elxy = sqrt(Edx^2+Edy^2);
ElxyScale = lTP/Elxy;
pts4 = pts;
pts4(2,1) = pts4(1,1) + ElxyScale*Edx;
pts4(2,2) = pts4(1,2) + ElxyScale*Edy;

%Find radius based on number of zeros before data at theta = 0
for i = 1:profWidth
        ptsTemp = pts4;
        ptsTemp(2,1) = pts4(1,1)+(pts4(2,1)-pts4(1,1))*cos((pi/180)*(i-profWidth/2))+(pts4(2,2)-pts4(1,2))*sin((pi/180)*(i-profWidth/2));
        ptsTemp(2,2) = pts4(1,2)-(pts4(2,1)-pts4(1,1))*sin((pi/180)*(i-profWidth/2))+(pts4(2,2)-pts4(1,2))*cos((pi/180)*(i-profWidth/2));
        clear profTemp
        profTemp = improfile(vq3,ptsTemp(:,1),ptsTemp(:,2),'bilinear');
    radI(i,1) = find(profTemp(:,1),1,'first');
    profTemp(profTemp==0,:) = []; %Removes Zeros to start profile right at shape edge
    dispProfXYThck(1:size(profTemp,1),i)= profTemp(:,1);
    plot(ptsTemp(1:2,1),ptsTemp(1:2,2))
    
end

if strcmp(type,'pill') == 1
    Prad = mean(radI);
end

%Set arc length at shape edge to sample
arcLength = 5; %microns
arcL = arcLength/.1625; %convert to pixels

%Radius of feature in microns
if strcmp(type,'circle') == 1
    rad = Crad;
elseif strcmp(type,'pill') == 1
    rad = Prad;
else
    rad = 0;
end

%Set the Number of Divisions Based on Arc Length
ArcTotal = rad*pi;
divisions = round(ArcTotal/arcL);
degreeTotal = ((arcL*divisions)/ArcTotal)*180;
divSize = degreeTotal/divisions;
dStart = degreeTotal/2;





%Pair up opposing slices (i.e. opposite angle relative to horizon)
for i = 1:round(divisions/2)
    pairs(i,1) = i;
    pairs(i,2) = (divisions+1)-i;
end

%Create colormap for associating profiles to pairs of slices
colormap = brewermap(size(pairs,1),'Dark2');

%% Divide shape radially and take the average profile for each division
% Done by averaging many profiles taken at 1degree increments within each
% division
clear Ys2
RadialProfsO = figure;
imshow(image.Trans,[]) % Show the heatmap for overlaying profiles
hold on
for j = 1:divisions
    [cidx, ~] = find(pairs == j, 1, 'first');
    clear dispProfXYThck
    for i = 1:round(divSize)
        
        %Positive direction trace
        ptsTemp = pts4;
        ptsTemp(2,1) = pts4(1,1)+(pts4(2,1)-pts4(1,1))*cos((pi/180)*(-dStart+i+(j-1)*divSize))+(pts4(2,2)-pts4(1,2))*sin((pi/180)*(-dStart+i+(j-1)*divSize));
        ptsTemp(2,2) = pts4(1,2)-(pts4(2,1)-pts4(1,1))*sin((pi/180)*(-dStart+i+(j-1)*divSize))+(pts4(2,2)-pts4(1,2))*cos((pi/180)*(-dStart+i+(j-1)*divSize));
        clear profTemp
        profTemp = improfile(vq3,ptsTemp(:,1),ptsTemp(:,2),'bilinear');
        profTemp(profTemp==0,:) = []; %Removes Zeros to start profile right at shape edge
        dispProfXYThck(1:size(profTemp,1),i)= profTemp(:,1);
        plot(ptsTemp(1:2,1),ptsTemp(1:2,2),'Color',colormap(cidx,:),'LineWidth',.5)
        
    end
    dispProfXYThckFinal = mean(dispProfXYThck,2,'omitnan');
    Ys2(1:size(dispProfXYThckFinal,1),j) = [dispProfXYThckFinal];
end
%%
%Save image of profile traces as overlay on displacement heatmap
title = ['\Profile Data\RadialProfileOverlay'];
savefile = [filePath title];
export_fig(RadialProfsO,savefile,'-native');

%Plot displacement as function of location on trace
Ys3 = Ys2;
Ys3(Ys3==0)=NaN;
for i = 1:size(pairs,1)
    Ys4(:,i) = mean([Ys3(:,pairs(i,1)),Ys3(:,pairs(i,2))],2,'omitnan');
end
Ys4(Ys4==0)=NaN;
Xs2 = 0:.1625:(.1625*(size(Ys4,1)-1));
RadialProfs = figure;
hold on
for i = 1:size(Ys4,2)
    [cidx, ~] = find(pairs == i, 1, 'first');
    plot(Xs2,Ys4(:,i),'Color',colormap(cidx,:))
end

%Save image of displacement profile plots
title = ['\Profile Data\RadialProfiles'];
savefile = [filePath title];
export_fig(RadialProfs,savefile,'-native');


end

