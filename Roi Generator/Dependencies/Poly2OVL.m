function output = Poly2OVL(filnm2,polygons)

%% Notes:
%write letters as uint8 (e.g. fwrite(file,uint8('A'),'uint8')
%%
filnm = '345Spline.rls';
pl=1;
fid=fopen(filnm);
%A(pl:pl+35) = fread(fid,36,'uint8=>char');
A(pl:pl+35,1) = fread(fid,36,'uint8');
pl=pl+36;
filnamnum2 = fread(fid,1,'uint32');
A(pl,1)=filnamnum2;
pl=pl+1;
A(pl:pl+filnamnum2,1)= fread(fid,filnamnum2+1,'uint8');
%A(pl:pl+filnamnum)= fread(fid,filnamnum+1,'uint8=>char')
pl=pl+filnamnum2+1;
numelements = fread(fid,1,'int32');
A(pl,1)=numelements;
pl=pl+1;
totsize = fread(fid,1,'int32');
A(pl,1)=totsize;
pl=pl+1;
A(pl:pl+20,1) = fread(fid,21,'int32');
pl=pl+21;
A(pl:pl+31,1) = fread(fid,32,'uint16');
%A(pl:pl+31) = fread(fid,32,'uint16=>char')
pl=pl+32;
A(pl:pl+9,1) = fread(fid,10,'int32');
pl=pl+10;

for i=1%:numelements
    %LOOP!
    %keyboard
    %for i=1:numelements
    drawelemtype = fread(fid,1,'int32');
    A(pl,1)=drawelemtype;
    pl=pl+1;
    %if drawelemtype~=18
    %    pl=pl-1;
    %    continue
    %end
    
    %Now I'm assuming it's polyline
    elemsize = fread(fid,1,'int32');
    A(pl,1)= elemsize;
    pl=pl+1;
    A(pl:pl+1,1) = fread(fid,2,'int32');
    pl=pl+2;
    A(pl:pl+1,1) = fread(fid,2,'float64');
    pl=pl+2;
    A(pl:pl+16,1) = fread(fid,17,'int32');
    pl=pl+17;
    A(pl:pl+31,1) = fread(fid,32,'uint16');
    %A(pl:pl+31) = fread(fid,32,'uint16=>char')
    pl=pl+32;
    A(pl,1) = fread(fid,1,'int32');
    pl=pl+1;
    A(pl:pl+8,1) = fread(fid,9,'int32');
    pl=pl+9;
    
    numvertex = fread(fid,1,'int32');
    A(pl,1)=numvertex;
    pl=pl+1;
    
    clear xxx
    clear yyy
    for j=1:numvertex
        xxx(i,j) = fread(fid,1,'float64');
        A(pl,1)=xxx(i,j);
        pl=pl+1;
        yyy(i,j) = fread(fid,1,'float64');
        A(pl,1)=yyy(i,j);
        pl=pl+1;
    end
end
A(pl,1) = fread(fid,1,'int32');
fclose(fid);

%%
%Total size of document (=64601+sum(elements)).
drawnumber = size(polygons,1);
elemsizes = zeros(drawnumber,2);
for i = 1:drawnumber
elemsizes(i,1) = size(polygons{i,1},1);
elemsizes(i,2) = 208+(16*elemsizes(i,1));
end
totalsize = 200+(sum(elemsizes(:,2)));

%%
pl=1;
fid=fopen(filnm2,'w');
fwrite(fid,A(pl:pl+35,1),'uint8');
pl=pl+36;
filnamnum = size(filnm2,2);
fwrite(fid,uint32(filnamnum),'uint32');
pl=pl+1;

for i = 1:size(filnm2,2)
    fwrite(fid,uint8(filnm2(i)),'uint8');
end

fwrite(fid,A(pl+filnamnum2,1),'uint8');
pl=pl+filnamnum2+1;
%Num Drawing Elements
fwrite(fid,int32(drawnumber),'int32');
pl=pl+1;
%Size
fwrite(fid,int32(totalsize),'int32');
pl=pl+1;
fwrite(fid,A(pl:pl+20,1),'int32');
pl=pl+21;
fwrite(fid,A(pl:pl+31,1),'uint16');
pl=pl+32;
fwrite(fid,A(pl:pl+9,1),'int32');
pl=pl+10;
plremember=pl;
for i=1:drawnumber
    numVertex = elemsizes(i,1);
    pl=plremember;
    fwrite(fid,int32(20),'int32');
    pl=pl+1;
    fwrite(fid,int32(elemsizes(i,2)),'int32');
    pl=pl+1;
    fwrite(fid,A(pl:pl+1,1),'int32');
    pl=pl+2;
    fwrite(fid,A(pl:pl+1,1),'float64');
    pl=pl+2;
    fwrite(fid,A(pl:pl+16,1),'int32');
    pl=pl+17;
    fwrite(fid,A(pl:pl+31,1),'uint16');
    pl=pl+32;
    fwrite(fid,A(pl,1),'int32');
    pl=pl+1;
    fwrite(fid,A(pl:pl+8,1),'int32');
    pl=pl+9;
    
    fwrite(fid,int32(numVertex),'int32');
    for j = 1:numVertex
    fwrite(fid,double(polygons{i,1}(j,1)),'float64');%x coordinate of vertex
    fwrite(fid,double(polygons{i,1}(j,2)),'float64');%y coordinate of vertex
    end
end
stest=size(A);
stest2=stest(1);
fwrite(fid,A(stest2,1),'int32');
fclose(fid);


output=0;