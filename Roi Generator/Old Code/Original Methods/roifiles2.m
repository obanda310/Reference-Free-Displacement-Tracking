function output = roifiles2(filnm,filnm2,custnum1,custnum2,custnum3,custnum4,drawnumber,splinex,spliney)

%%
filnm = '345Spline.rls';

pl=1;

fid=fopen(filnm);
%A(pl:pl+35) = fread(fid,36,'uint8=>char');
A(pl:pl+35,1) = fread(fid,36,'uint8');
pl=pl+36;
filnamnum = fread(fid,1,'uint32');
A(pl,1)=filnamnum;
pl=pl+1;
A(pl:pl+filnamnum,1)= fread(fid,filnamnum+1,'uint8');
%A(pl:pl+filnamnum)= fread(fid,filnamnum+1,'uint8=>char')
pl=pl+filnamnum+1;
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
%%
A(pl,1) = fread(fid,1,'int32');

fclose(fid);
%%
pl=1;
fid=fopen(filnm2,'w');
fwrite(fid,A(pl:pl+35,1),'uint8');
pl=pl+36;
fwrite(fid,uint32(7),'uint32');
pl=pl+1;
fwrite(fid,uint8('J'),'uint8');
fwrite(fid,uint8('C'),'uint8');
fwrite(fid,uint8('C'),'uint8');
fwrite(fid,uint8(custnum1),'uint8');
fwrite(fid,uint8(custnum2),'uint8');
fwrite(fid,uint8(custnum3),'uint8');
fwrite(fid,uint8(custnum4),'uint8');

%fwrite(fid,uint8('C'),'uint8');
%fwrite(fid,uint8('u'),'uint8');
%fwrite(fid,uint8('s'),'uint8');
%fwrite(fid,uint8('t'),'uint8');
%fwrite(fid,uint8('o'),'uint8');
%fwrite(fid,uint8('m'),'uint8');
%fwrite(fid,uint8(custnum1),'uint8');



fwrite(fid,A(pl+filnamnum,1),'uint8');
pl=pl+filnamnum+1;
%Num Drawing Elements
fwrite(fid,int32(drawnumber),'int32');
pl=pl+1;
%Size
fwrite(fid,int32((200+(drawnumber*272))),'int32');
pl=pl+1;
fwrite(fid,A(pl:pl+20,1),'int32');
pl=pl+21;
fwrite(fid,A(pl:pl+31,1),'uint16');
pl=pl+32;
fwrite(fid,A(pl:pl+9,1),'int32');
pl=pl+10;
plremember=pl;
for i=1:drawnumber
    pl=plremember;
    fwrite(fid,int32(20),'int32');
    pl=pl+1;
    fwrite(fid,int32(272),'int32');
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
    
    fwrite(fid,int32(4),'int32');
    fwrite(fid,double(splinex(i,1)),'float64');
    fwrite(fid,double(spliney(i,1)),'float64');
    fwrite(fid,double(splinex(i,2)),'float64');
    fwrite(fid,double(spliney(i,2)),'float64');
    fwrite(fid,double(splinex(i,3)),'float64');
    fwrite(fid,double(spliney(i,3)),'float64');
    fwrite(fid,double(splinex(i,4)),'float64');
    fwrite(fid,double(spliney(i,4)),'float64');
end
stest=size(A);
stest2=stest(1);
fwrite(fid,A(stest2,1),'int32');
fclose(fid);


output=0;