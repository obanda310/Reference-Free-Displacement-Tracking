function output = Poly2Regions(filnm2,polygons)

%%
filnm = 'Jack15.Regions';
load('regionIndex.mat')
fid = fopen(filnm);
fid2 = fopen(filnm);
Alength = size(fread(fid2));

%%  Read the Header
h1=1;
H(h1:h1+89,1) = fread(fid,90,'uint8');
h1=h1+90;
totsize = fread(fid,1,'int32');
H(h1,1)=totsize;
h1=h1+1;
H(h1:h1+277,1)= fread(fid,278,'uint8');
h1=h1+278;
for i = 2:size(regionsIndcs)-1
        H(h1,1) = fread(fid,1,'int32');
        h1 = h1+1;
        if regionsIndcs(i+1,1)-regionsIndcs(i,1)~=4
            H(h1:h1+(regionsIndcs(i+1,1)-regionsIndcs(i,1))-5,1) = fread(fid,(regionsIndcs(i+1,1)-regionsIndcs(i,1))-4,'uint8');
            h1=h1+(regionsIndcs(i+1,1)-regionsIndcs(i,1))-4;
        end
end
H(h1,1) = fread(fid,1,'int32');
h1=h1+1;
H(h1:h1+(1362-1178)-1,1) = fread(fid,1362-1178,'uint8');
h1 = h1+(1362-1178);

%% Read the Elements
clear E
e1 = 1;
%size(fread(fid,'uint8'))
    %LOOP!
    %keyboard
    %for i=1:numelements
    drawelemtype = fread(fid,1,'int32');
    E(e1,1)=drawelemtype;
    e1=e1+1;
    %if drawelemtype~=18
    %    p1=p1-1;
    %    continue
    %end
    
    %Now I'm assuming it's polyline
    elemsize = fread(fid,1,'int32');
    E(e1,1)= elemsize;
    e1=e1+1;
    E(e1:e1+1,1) = fread(fid,2,'int32');
    e1=e1+2;
    E(e1:e1+1,1) = fread(fid,2,'float64');
    e1=e1+2;
    E(e1:e1+16,1) = fread(fid,17,'int32');
    e1=e1+17;
    E(e1:e1+31,1) = fread(fid,32,'uint16');
    %A(p1:p1+31) = fread(fid,32,'uint16=>char')
    e1=e1+32;
    E(e1,1) = fread(fid,1,'int32');
    e1=e1+1;
    E(e1:e1+8,1) = fread(fid,9,'int32');
    e1=e1+9;
    
    numvertex = fread(fid,1,'int32');
    E(e1,1)=numvertex;
    e1=e1+1;
    
    clear xxx
    clear yyy
    for j=1:numvertex
        xxx(i,j) = fread(fid,1,'float64');
        E(e1,1)=xxx(i,j);
        e1=e1+1;
        yyy(i,j) = fread(fid,1,'float64');
        E(e1,1)=yyy(i,j);
        e1=e1+1;
    end
    
%% Read the Footer
clear F f1

f1=1;
numelements = 2286;
skip = (numelements-1)*272;
Askip = fread(fid,skip,'uint8');

F(f1:f1+63238,1) = fread(fid,63239,'uint8');
fclose(fid);

%%

%%Use fragments in 'H, E, F' to generate a .Regions file
h1=1;
fid=fopen(filnm2,'w');
%Total size of document (=64601+sum(elements)).
drawnumber = size(polygons,1);
elemsizes = zeros(drawnumber,2);
for i = 1:drawnumber
elemsizes(i,1) = size(polygons{i,1},1);
elemsizes(i,2) = 208+(16*elemsizes(i,1));
end
totalsize = 64601+(sum(elemsizes(:,2)));
disp(['Regions file has ' num2str(sum(elemsizes(:,1))) ' vertices'])
%% Write the Header

fwrite(fid,H(h1:h1+89,1),'uint8');
h1=h1+90;
fwrite(fid,int32(totalsize),'int32');
h1=h1+1;
fwrite(fid,H(h1:h1+277,1),'uint8');
h1=h1+278;

% There are several index values here in the header that point to specific
% locations in the footer. To be recalculated, do totalsize-index. indices
% should be saved in the regionIndex.mat dependency file.

for i = 2:22
        fwrite(fid,int32(totalsize-regionsIndcs(i,2)-1),'int32');
        h1 = h1+1;
        if regionsIndcs(i+1,1)-regionsIndcs(i,1)~=4
            fwrite(fid,H(h1:h1+(regionsIndcs(i+1,1)-regionsIndcs(i,1))-5,1),'uint8');
            h1=h1+(regionsIndcs(i+1,1)-regionsIndcs(i,1))-4;
        end
end

fwrite(fid,int32(totalsize-4),'int32');
h1=h1+1;
fwrite(fid,int32(totalsize-3),'int32');
h1=h1+1;
fwrite(fid,int32(totalsize-2),'int32');
h1=h1+1;
fwrite(fid,int32(1),'int32');
h1=h1+1;
fwrite(fid,H(h1:h1+(regionsIndcs(27,1)-regionsIndcs(26,1))-5,1),'uint8');
h1=h1+(regionsIndcs(27,1)-regionsIndcs(26,1));

%Number of Drawing Elements
fwrite(fid,int32(drawnumber),'int32');
h1=h1+1;
fwrite(fid,int32(totalsize-regionsIndcs(28,2)),'int32');
h1=h1+1;
fwrite(fid,int32(3),'int32');
h1=h1+1;
fwrite(fid,int32(0),'int32'); %lost some zeros for some reason. This fixes
fwrite(fid,int8(0),'int8');
h1=h1+1;
fwrite(fid,H(h1:h1+(1357-1178)-1,1),'uint8');
%h1 = h1+(1357-1178);
%fwrite(fid,int16(0),'uint16');

%% Write the Elements
for i=1:drawnumber
    numVertex = elemsizes(i,1);
    e1=1;
    fwrite(fid,int32(20),'int32'); % "20" is polyline I think
    e1=e1+1;
    fwrite(fid,int32(elemsizes(i,2)),'int32');%Number of bits per element(=(numVertex*16)+208)
    e1=e1+1;
    fwrite(fid,E(e1:e1+1,1),'int32');
    e1=e1+2;
    fwrite(fid,E(e1:e1+1,1),'float64');
    e1=e1+2;
    fwrite(fid,E(e1:e1+16,1),'int32');
    e1=e1+17;
    fwrite(fid,E(e1:e1+31,1),'uint16');
    e1=e1+32;
    fwrite(fid,E(e1,1),'int32');
    e1=e1+1;
    fwrite(fid,E(e1:e1+8,1),'int32');
    e1=e1+9;
    
    fwrite(fid,int32(numVertex),'int32');
    for j = 1:numVertex
    fwrite(fid,double(polygons{i,1}(j,1)),'float64');%x coordinate of vertex
    fwrite(fid,double(polygons{i,1}(j,2)),'float64');%y coordinate of vertex
    end
end

%% Write the Footer
f1=1;
fwrite(fid,F(f1+9:end,1),'uint8');
fclose(fid);

output=0;