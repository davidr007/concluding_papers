%function B = sacmodpo(filnm,filnm2,entypcurrent,t0,t1)
%modifies according to selected time window t0 to t1
%t0 and t1 are points NOT SECONDS
%Erdinc Saygin 
%july 2004

function B = sacmodpo(filnm,filnm2,entypcurrent,t0,t1)


[fid,message] = fopen(filnm,'r',entypcurrent);
if fid==-1
    display('could not find the file');
    return;
end






[farr,count] = fread(fid,70,'float32');
[iarr,count2] = fread(fid,40,'int32');
[karr,count3] = fread(fid,192,'schar');

karr=karr';





%header values
 delta=farr(1);
% depmin=farr(2);
% depmax=farr(3);
% scale=farr(4);
% odelta=farr(5);
 b=farr(6);
 e=farr(7);
% o=farr(8);
% a=farr(9);
% t0=farr(11);
% t1=farr(12);
% t2=farr(13);
% t3=farr(14);
% t4=farr(15);
% t5=farr(16);
% t6=farr(17);
% t7=farr(18);
% t8=farr(19);
% t9=farr(20);
% f=farr(21);
% resp0=farr(22);
% resp1=farr(23);
% resp2=farr(24);
% resp3=farr(25);
% resp4=farr(26);
% resp5=farr(27);
% resp6=farr(28);
% resp7=farr(29);
% resp8=farr(30);
% resp9=farr(31);
% stla=farr(32);
% stlo=farr(33);
% stel=farr(34);
% stdp=farr(35);
% evla=farr(36);
% evlo=farr(37);
% evel=farr(38);
% evdp=farr(39);
% mag=farr(40);
% user0=farr(41);
% user1=farr(42);
% user2=farr(43);
% user3=farr(44);
% user4=farr(45);
% user5=farr(46);
% user6=farr(47);
% user7=farr(48);
% user8=farr(49);
% user9=farr(50);
% dist=farr(51);
% az=farr(52);
% baz=farr(53);
% gcarc=farr(54);
% depmen=farr(55);
% cmpaz=farr(56);
% cmpinc=farr(57);
% xminimum=farr(58);
% xmaximum=farr(59);
% yminimum=farr(60);
% ymaximum=farr(61);
% %integer values
% nzyear=iarr(1);
% nzjday=iarr(2);
% nzhour=iarr(3);
% nzmin=iarr(4);
% nzsec=iarr(5);
% nzmsec=iarr(6);
% nvhdr=iarr(7);
% norid=iarr(8);
% nevid=iarr(9);
npts=iarr(10);
% nwfid=iarr(12);
% nxsize=iarr(13);
% nysize=iarr(14);
% iftype=iarr(16);
% idep=iarr(17);
% iztype=iarr(18);
% iinst=iarr(20);
% istreg=iarr(21);
% ievreg=iarr(22);
% ievtyp=iarr(23);
% iqual=iarr(24);
% isynth=iarr(25);
% imagtyp=iarr(26);
% imagsrc=iarr(27);
% %logical part
% leven=iarr(36);
% lpspol=iarr(37);
% lovrok=iarr(38);
% lcalda=iarr(39);
% %-----------
% kstnm=char(karr(1:8));
% kevnm=char(karr(9:24));
% khole=char(karr(25:32));
% ko=char(karr(33:40));
% ka=char(karr(41:48));
% kt0=char(karr(49:56));
% kt1=char(karr(57:64));
% kt2=char(karr(65:72));
% kt3=char(karr(73:80));
% kt4=char(karr(81:88));
% kt5=char(karr(89:96));
% kt6=char(karr(97:104));
% kt7=char(karr(105:112));
% kt8=char(karr(113:120));
% kt9=char(karr(121:128));
% kf=char(karr(129:136));
% kuser0=char(karr(137:144));
% kuser1=char(karr(145:152));
% kuser2=char(karr(153:160));
% kcmnp=char(karr(161:168));
% knetwk=char(karr(169:176));
% kdatrd=char(karr(177:184));
% kinst=char(karr(185:192));
%-----reading data
data=fread(fid,npts,'float32');

%-----extracting the wanted part from main data array

beg=t0; %fix((t0-b)./delta)
endd=t1; %((t1-b)./delta)

npts2=endd-beg+1;
data2=data(beg:endd);
%----updating the necessary header values
% farr(6)=b+(t0.*delta);%b
% farr(7)=b+(t1.*delta); %e
farr(6)=b+((t0-1).*delta);%b
farr(7)=b+((t1-1).*delta); %e



iarr(10)=npts2;




%filnm2=strcat(filnm,'_mod');
[fid2,message] = fopen(filnm2,'w',entypcurrent);
fwrite(fid2,farr,'float32');
fwrite(fid2,iarr,'int32');
fwrite(fid2,karr,'schar');
fwrite(fid2,data2,'float32');

fclose('all');
