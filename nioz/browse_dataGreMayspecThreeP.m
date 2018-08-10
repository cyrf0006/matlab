clear
close all

%%
load /media/Seagate1TB/NIOZ/thermistordata/GRE12may/GreMay2012_1.txt
ctd=GreMay2012_1;
[nctd,mctd]=size(ctd);

%%
cruise_path='/media/Seagate1TB/NIOZ/thermistordata/GRE12may/';
thermistor_path=['thermistors'];

calib_path=fullfile('/media/Seagate1TB/NIOZ/calibrationdata');
calib_name='S_ijkbad_2011_03_14_range15_22';
load(fullfile(calib_path,calib_name));

wrong=[];true=[];

t=thermistor();
strnumber=num2str(reindex(101,wrong,true),'%04d');
filename=fullfile(cruise_path,thermistor_path,[strnumber '_000.SEN.dec.bin']);
t=append_data(t,filename,'nodata');
t_for_range=ref_date(t,datenum(2001,01,01));

dag1=149.4219;%in water...10...degr3
dag2=152.3018;
% dag1=149.5;%in water...10...degr3
% dag2=150.5;

% - FC (useless)
% $$$ ra=[8 19];
% $$$ cont=20;%...1/cont degr

range_wanted=[dag1 dag2];
range=range_estimator(t_for_range,range_wanted,'yeardays')
skipcal=50;
skipdata=5;

%part 1
thres=0.01;
skip=skipcal;
s=string();
for p=[101:102 104:109]
  clear t
  strnumber=num2str(reindex(p,wrong,true),'%04d');
  filename=fullfile(cruise_path,thermistor_path,[strnumber '_000.SEN.dec.bin']);
  t=thermistor();
      if exist(filename)
          t=append_data(t,filename,'skip',skip,'range',range,'noaccel');
      end
  s=append_therm(s,t);
end
%%
s=ref_date(s,datenum(2001,1,1));
eval(['s=import_calib(s,' calib_name ');']);
calib_path=fullfile('/media/Seagate1TB/NIOZ/calibrationdata');
calib_name='S_ijkbad_2010_12_07';
load(fullfile(calib_path,calib_name));
s=ref_date(s,datenum(2001,1,1));
eval(['s=import_calib(s,' calib_name ');']);

s=auto_scale(s,'treshold',thres);%bathcalib
s_calib=remove_data(s);

skip=skipdata;
s=string();
for p=[101:102 104:109]%[invoegen 69:62]
  clear t
  strnumber=num2str(reindex(p,wrong,true),'%04d');
  filename=fullfile(cruise_path,thermistor_path,[strnumber '_000.SEN.dec.bin']);
  t=thermistor();
      if exist(filename)
          t=append_data(t,filename,'skip',skip,'range',range,'noaccel');
      end
  s=append_therm(s,t);
end
%%
s=ref_date(s,datenum(2001,1,1));
s=import_calib(s,s_calib);

[tijd,TT]=string2mat(s,'calib','yearday');
T1=TT;
[mt1,nt]=size(T1);

%part 2
% t=thermistor();
% strnumber=num2str(reindex(101,wrong,true),'%04d');
% filename=fullfile(cruise_path,thermistor_path,[strnumber '_000.SEN.dec.bin']);
% t=append_data(t,filename,'nodata');
% t_for_range=ref_date(t,datenum(2001,01,01));
% range=range_estimator(t_for_range,range_wanted,'yeardays')
% clear s
% clear s_calib

thres=0.01;
skip=skipcal;
s=string();
for p=[110:147]
  clear t
  strnumber=num2str(reindex(p,wrong,true),'%04d');
  filename=fullfile(cruise_path,thermistor_path,[strnumber '_000.SEN.dec.bin']);
  t=thermistor();
      if exist(filename)
          t=append_data(t,filename,'skip',skip,'range',range,'noaccel');
      end
  s=append_therm(s,t);
end
%%
%calib_path=fullfile('D:\NIOZHST\bathcalib\ijkbad_14_03_2011\data');
calib_path=fullfile('/media/Seagate1TB/NIOZ/calibrationdata');
calib_name='S_ijkbad_2011_03_14_range15_22';
load(fullfile(calib_path,calib_name));
s=ref_date(s,datenum(2001,1,1));
eval(['s=import_calib(s,' calib_name ');']);
calib_path=fullfile('/media/Seagate1TB/NIOZ/calibrationdata');
%calib_path=fullfile('D:\NIOZHST\bathcalib\ijkbad_07_12_2010\data');
calib_name='S_ijkbad_2010_12_07';
load(fullfile(calib_path,calib_name));
s=ref_date(s,datenum(2001,1,1));
eval(['s=import_calib(s,' calib_name ');']);

s=auto_scale(s,'treshold',thres);%bathcalib
s_calib=remove_data(s);

skip=skipdata;
s=string();
for p=[110:147]
  clear t
  strnumber=num2str(reindex(p,wrong,true),'%04d');
  filename=fullfile(cruise_path,thermistor_path,[strnumber '_000.SEN.dec.bin']);
  t=thermistor();
      if exist(filename)
          t=append_data(t,filename,'skip',skip,'range',range,'noaccel');
      end
  s=append_therm(s,t);
end
%%
s=ref_date(s,datenum(2001,1,1));
s=import_calib(s,s_calib);

[tijd,TT]=string2mat(s,'calib','yearday');

T=TT;

for p=[19]
        TT(p,:)=mean(TT(p-1:2:p+1,:));
end
T2=TT;
[mt2,nt]=size(T2);

%part 3

thres=0.001;
skip=skipcal;
s=string();
for p=[148:175 177:190 192:195]
  clear t
  strnumber=num2str(reindex(p,wrong,true),'%04d');
  filename=fullfile(cruise_path,thermistor_path,[strnumber '_000.SEN.dec.bin']);
  t=thermistor();
      if exist(filename)
          t=append_data(t,filename,'skip',skip,'range',range,'noaccel');
      end
  s=append_therm(s,t);
end
%%
calib_path=fullfile('/media/Seagate1TB/NIOZ/calibrationdata');
%calib_path=fullfile('D:\NIOZHST\bathcalib\ijkbad_14_03_2011\data');
calib_name='S_ijkbad_2011_03_14_range15_22';
load(fullfile(calib_path,calib_name));
s=ref_date(s,datenum(2001,1,1));
eval(['s=import_calib(s,' calib_name ');']);
calib_path=fullfile('/media/Seagate1TB/NIOZ/calibrationdata');
%calib_path=fullfile('D:\NIOZHST\bathcalib\ijkbad_07_12_2010\data');
calib_name='S_ijkbad_2010_12_07';
load(fullfile(calib_path,calib_name));
s=ref_date(s,datenum(2001,1,1));
eval(['s=import_calib(s,' calib_name ');']);

s=auto_scale(s,'treshold',thres);%bathcalib
s_calib=remove_data(s);

skip=skipdata;
s=string();
for p=[148:175 177:190 192 001 193:195]%[invoegen 69:62]
  clear t
  strnumber=num2str(reindex(p,wrong,true),'%04d');
  filename=fullfile(cruise_path,thermistor_path,[strnumber '_000.SEN.dec.bin']);
  t=thermistor();
      if exist(filename)
          t=append_data(t,filename,'skip',skip,'range',range,'noaccel');
      end
  s=append_therm(s,t);
end
%%
s=ref_date(s,datenum(2001,1,1));
s=import_calib(s,s_calib);

[tijd,TT]=string2mat(s,'calib','yearday');
T=TT;
for p=[18 38]
        TT(p,:)=mean(TT(p-1:2:p+1,:));
end
for p=[44]
        TT(p,:)=TT(p-1,:)*2/3+TT(p+2,:)/3;
        TT(p+1,:)=TT(p-1,:)/3+TT(p+2,:)*2/3;
end
T3=TT;
[mt3,nt]=size(T3);

TT=zeros(mt1+mt2+mt3,nt);
TT(1:mt1,:)=T1;
TT(mt1+1:mt1+mt2,:)=T2;
TT(mt1+mt2+1:mt1+mt2+mt3,:)=T3;
T=flipdim(TT,1);
[mt,nt]=size(T);
dz=-[0:mt-1]*0.333+33.2;
nf1=1;
nf2=nt;
plot(dz,mean(T'),ctd(:,1),ctd(:,3),'g',ctd(:,1)-3,ctd(:,3),'g--')
% break

drdT=0.4188*1e-2;
deltaz=0.3333;
  [Tsort,I]=sort(T,1);
  imagesc(tijd,-dz,TT)
  break
%%
%spectra

% tap=flipdim(Tsort,1);
i01=TT(2,:);
w01=TT(40,:);
v01=TT(70,:);
u01=TT(70,:);
v02=TT(90,:);
u02=TT(90,:);
m1=spectr(i01,i01,40000,20000);
m2=spectr(w01,w01,40000,20000);
m3=spectr(u01,v01,40000,20000);
m4=spectr(u02,v02,40000,20000);
m1=spectr(i01,i01,4000,2000);
m2=spectr(w01,w01,4000,2000);
m3=spectr(u01,v01,4000,2000);
m4=spectr(u02,v02,4000,2000);

figure(5)
clf
fac=55;
B='';
seIRZoEI(m3,m4,m2/10,m1,fac,86400/(skipdata+1),B)
print -djpeg Gre12Tspec

[mm,nn]=size(T);
% maT=T;
% mst=5;
% for ii=1+mst:nn-mst
%     for jj=1:mm
%         maT(jj,ii)=mean(T(jj,ii-mst:ii+mst));
%     end
% end

vw=12;
ddz=-dz;
% vcTs=zeros(nn,1);
% for ii=1:nn
%     ff=find(Tsort(:,ii)<vw(1));
%     [k1,k2]=size(ff);
%     vcTs(ii)=ddz(k1+1)+(vw(1)-Tsort(k1,ii))/(Tsort(k1+1,ii)-Tsort(k1,ii))*deltaz;
% end
% vcT=zeros(nn,1);
% ddz=-dz;
% for ii=1:nn
%     ff=find(maT(:,ii)<vw(1));
%     [k1,k2]=size(ff);
%     vcT(ii)=ddz(k1+1)+(vw(1)-maT(k1,ii))/(maT(k1+1,ii)-maT(k1,ii))*deltaz;
% end
vcT=zeros(nn,1);
for ii=1:nn
    ff=find(T(:,ii)<vw(1));
    [k1,k2]=size(ff);
    vcT(ii)=ddz(k1+1)+(vw(1)-T(k1,ii))/(T(k1+1,ii)-T(k1,ii))*deltaz;
end
gvcT=2/(skipdata+1)*gradient(vcT);
i01=gvcT(23:4.9e4);
w01=TT(40,:);
v01=TT(70,:);
u01=TT(70,:);
v02=TT(90,:);
u02=TT(90,:);
% m1=spectr(i01,i01,40000,20000);
% m2=spectr(w01,w01,40000,20000);
% m3=spectr(u01,v01,40000,20000);
% m4=spectr(u02,v02,40000,20000);
m1=spectr(i01,i01,4000,2000);
m2=spectr(w01,w01,4000,2000);
m3=spectr(u01,v01,4000,2000);
m4=spectr(u02,v02,4000,2000);

figure(5)
clf
fac=55;
B='';
seIRZoEI(m3,m4,m2,m1,fac,2*86400/(skipdata+1),B)
% print -djpeg Gre12TspecfromW

