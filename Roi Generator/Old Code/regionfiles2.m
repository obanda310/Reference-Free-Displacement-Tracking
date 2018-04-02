function output = regionfiles2(filnm,filnm2,custnum1,custnum2,custnum3,custnum4,drawnumber,sp1inex,sp1iney)

%%
%keyboard
clear A
clear ver1
clear hor1
clear ver2
clear hor2
load('regionIndex.mat')
p1=1;

fid = fopen(filnm);
fid2 = fopen(filnm);
Alength = size(fread(fid2));


A(p1:p1+89,1) = fread(fid,90,'uint8');
p1=p1+90;
totsize = fread(fid,1,'int32');
A(p1,1)=totsize;
p1=p1+1;
A(p1:p1+277,1)= fread(fid,278,'uint8');
p1=p1+278;
for i = 2:size(regionsIndcs)-1
        A(p1,1) = fread(fid,1,'int32');
        p1 = p1+1;
        if regionsIndcs(i+1,1)-regionsIndcs(i,1)~=4
            A(p1:p1+(regionsIndcs(i+1,1)-regionsIndcs(i,1))-5,1) = fread(fid,(regionsIndcs(i+1,1)-regionsIndcs(i,1))-4,'uint8');
            p1=p1+(regionsIndcs(i+1,1)-regionsIndcs(i,1))-4;
        end
end
A(p1,1) = fread(fid,1,'int32');
p1=p1+1;
A(p1:p1+(1362-1178)-1,1) = fread(fid,1362-1178,'uint8');
p1 = p1+(1362-1178);
%%
%size(fread(fid,'uint8'))
    %LOOP!
    %keyboard
    %for i=1:numelements
    drawelemtype = fread(fid,1,'int32');
    A(p1,1)=drawelemtype;
    p1=p1+1;
    %if drawelemtype~=18
    %    p1=p1-1;
    %    continue
    %end
    
    %Now I'm assuming it's polyline
    elemsize = fread(fid,1,'int32');
    A(p1,1)= elemsize;
    p1=p1+1;
    A(p1:p1+1,1) = fread(fid,2,'int32');
    p1=p1+2;
    A(p1:p1+1,1) = fread(fid,2,'float64');
    p1=p1+2;
    A(p1:p1+16,1) = fread(fid,17,'int32');
    p1=p1+17;
    A(p1:p1+31,1) = fread(fid,32,'uint16');
    %A(p1:p1+31) = fread(fid,32,'uint16=>char')
    p1=p1+32;
    A(p1,1) = fread(fid,1,'int32');
    p1=p1+1;
    A(p1:p1+8,1) = fread(fid,9,'int32');
    p1=p1+9;
    
    numvertex = fread(fid,1,'int32');
    A(p1,1)=numvertex;
    p1=p1+1;
    
    clear xxx
    clear yyy
    for j=1:numvertex
        xxx(i,j) = fread(fid,1,'float64');
        A(p1,1)=xxx(i,j);
        p1=p1+1;
        yyy(i,j) = fread(fid,1,'float64');
        A(p1,1)=yyy(i,j);
        p1=p1+1;
    end

numelements = 2286;
skip = (numelements-1)*272;
Askip = fread(fid,skip,'uint8');


A(p1:p1+63238,1) = fread(fid,63239,'uint8');
fclose(fid);
%%

%%Use fragments in 'A' to generate a .Regions file
p1=1;
fid=fopen(filnm2,'w');
totalsize = 64601+(272*drawnumber);

%For Reference

% A(p1:p1+277,1)= fread(fid,278,'uint8');
% p1=p1+278;
% for i = 2:size(regionsIndcs)
%         A(p1,1) = fread(fid,1,'uint32');
%         p1 = p1+1;
%         if regionsIndcs(i+1,1)-regionsIndcs(i,1)~=4
%             A(p1:p1+(regionsIndcs(i+1,1)-regionsIndcs(i,1))-4,1) = fread(fid,(regionsIndcs(i+1,1)-regionsIndcs(i,1))-4,'uint8');
%             p1=p1+(regionsIndcs(i+1,1)-regionsIndcs(i,1))-4;
%         end
% end
% 
% A(p1:p1+(1362-1178)-1,1) = fread(fid,1362-1178,'uint8');
% p1 = p1+(1362-1178);

fwrite(fid,A(p1:p1+89,1),'uint8');
p1=p1+90;
fwrite(fid,int32(totalsize),'int32');
p1=p1+1;
fwrite(fid,A(p1:p1+277,1),'uint8');
p1=p1+278;
for i = 2:22
        fwrite(fid,int32(totalsize-regionsIndcs(i,2)-1),'int32');
        p1 = p1+1;
        if regionsIndcs(i+1,1)-regionsIndcs(i,1)~=4
            fwrite(fid,A(p1:p1+(regionsIndcs(i+1,1)-regionsIndcs(i,1))-5,1),'uint8');
            p1=p1+(regionsIndcs(i+1,1)-regionsIndcs(i,1))-4;
        end
end

fwrite(fid,int32(totalsize-4),'int32');
p1=p1+1;
fwrite(fid,int32(totalsize-3),'int32');
p1=p1+1;
fwrite(fid,int32(totalsize-2),'int32');
p1=p1+1;
fwrite(fid,int32(1),'int32');
p1=p1+1;
fwrite(fid,A(p1:p1+(regionsIndcs(27,1)-regionsIndcs(26,1))-5,1),'uint8');
p1=p1+(regionsIndcs(27,1)-regionsIndcs(26,1));

%Num Drawing Elements
fwrite(fid,int32(drawnumber),'int32');
p1=p1+1;
fwrite(fid,int32(totalsize-regionsIndcs(28,2)),'int32');
p1=p1+1;
fwrite(fid,int32(3),'int32');
p1=p1+1;
fwrite(fid,int32(0),'int32'); %lost some zeros for some reason. This fixes
fwrite(fid,int8(0),'int8');
p1=p1+1;
fwrite(fid,A(p1:p1+(1357-1178)-1,1),'uint8');
p1 = p1+(1357-1178);
%fwrite(fid,int16(0),'uint16');
p1remember=p1;
for i=1:drawnumber
    p1=p1remember;
    fwrite(fid,int32(20),'int32');
    p1=p1+1;
    fwrite(fid,int32(272),'int32');
    p1=p1+1;
    fwrite(fid,A(p1:p1+1,1),'int32');
    p1=p1+2;
    fwrite(fid,A(p1:p1+1,1),'float64');
    p1=p1+2;
    fwrite(fid,A(p1:p1+16,1),'int32');
    p1=p1+17;
    fwrite(fid,A(p1:p1+31,1),'uint16');
    p1=p1+32;
    fwrite(fid,A(p1,1),'int32');
    p1=p1+1;
    fwrite(fid,A(p1:p1+8,1),'int32');
    p1=p1+9;
    
    fwrite(fid,int32(4),'int32');
    fwrite(fid,double(sp1inex(i,1)),'float64');
    fwrite(fid,double(sp1iney(i,1)),'float64');
    fwrite(fid,double(sp1inex(i,2)),'float64');
    fwrite(fid,double(sp1iney(i,2)),'float64');
    fwrite(fid,double(sp1inex(i,3)),'float64');
    fwrite(fid,double(sp1iney(i,3)),'float64');
    fwrite(fid,double(sp1inex(i,4)),'float64');
    fwrite(fid,double(sp1iney(i,4)),'float64');
end
fwrite(fid,A(p1+9:end,1),'uint8');
% stest=size(A);
% stest2=stest(1);
% fwrite(fid,A(stest2,1),'int32');
fclose(fid);


output=0;