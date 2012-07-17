
%function A = sacwrite(filnm,entypcurrent,data,delta,npts,b,e,nzyear,nzjday,nzhour,nzmin,nzsec,kstnm,kcmnp,stla,stlo,evla,evlo,cmpaz,cmpinc)
%filnm=File name of data
%entypcurrent= current file type ('b' or 'l')

%'b' is for UNIX files (big endian)
%'l' is for PC files (little endian)
%Erdinc Saygin 
%july 2004

function A = sacwrite(filnm,entypcurrent,data,delta,npts,b,e,nzyear,nzjday,nzhour,nzmin,nzsec,kstnm,kcmnp,stla,stlo,evla,evlo,cmpaz,cmpinc)

[fid,message] = fopen(filnm,'wb',entypcurrent);

if fid==-1
    display('there is a problem, file could not be formed');
    return;
end

%initialize
farr(1:70)=-12345.0;
iarr(1:40)=-12345;

karr(1,1:192)=37;
karr=karr';

%-------------

farr(1)=delta;
farr(6)=b;
farr(7)=e;

farr(32)=stla;
farr(33)=stlo;
farr(36)=evla;
farr(37)=evlo;
farr(58)=cmpaz;
farr(59)=cmpinc;

iarr(1)=nzyear;
iarr(2)=nzjday;
iarr(3)=nzhour;
iarr(4)=nzmin;
iarr(5)=nzsec;

iarr(7)=6; %nvhdr;
iarr(10)=npts;
iarr(16)=1; %iftype=;


[sr sc]=size(kstnm);
[sr2 sc2]=size(kcmnp);
if(sc<8)
    kstnm(1,sc+1:8)=' ';
end
if(sc2<8)
    kcmnp(1,sc2+1:8)=' ';
end





karr(1:8)=uint8(kstnm(1,1:8));
karr(161:168)=uint8(kcmnp(1,1:8));

A=karr;


fwrite(fid,farr,'float32');
fwrite(fid,iarr,'int32'); %int
count=fwrite(fid,karr,'int8'); %schar di!!!
count2=fwrite(fid,data,'float32');

fclose(fid);

