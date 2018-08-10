%
%diepe minima zoeken (plot(tijd,Tsort(1:10,:))

%%%%
%NOTE!!
% turbul is goed, maar imagesc --> pcolor omdat vertikale schaal varieert
% 0.6 -- 1.0 m
%%%%%%

clear
close all
%%
load /home/cyrf0006/research/NIOZ/RockallBank/CTD360/ascii1m/64PE360_S93C01fil_ctm_drv_bav_buo.asc
aa=X64PE360_S93C01fil_ctm_drv_bav_buo;
figure;plot(aa(:,2),-aa(:,1))
ff=find(aa(:,1)>500);
[n,m]=size(ff);
k1=ff(1);
k2=ff(n);
figure;plot(aa(k1:k2,2),-aa(k1:k2,1))
ft1=polyfit(aa(k1:k2,2),aa(k1:k2,11),1)
ag1=polyval(ft1,aa(k1:k2,2));
figure;plot(aa(k1:k2,11),aa(k1:k2,2),'k.',ag1,aa(k1:k2,2),'r')
print -djpeg ROC12Tsig1000
%%
load /media/Seagate1TB/NIOZ/thermistordata/ROC12/roc12.mat
[nA,mA]=size(SerEA1cnt);
I1=SerEA1cnt*0.45;dI1=I1;mI1=mean(I1);
for ii=1:mA;dI1(:,ii)=I1(:,ii)-mI1(ii);end
pA=AnDepthmm/1000;
ww=SerVmmpersec/1000;
ya=273+SerDay+SerHour/24+SerMin/1440+SerSec/86400+SerHund/86400/100;
clear Ser*
clear An*
%%
cruise_path='/media/Seagate1TB/NIOZ/thermistordata/ROC12/';
thermistor_path=[''];

dummy=xlsread(fullfile(cruise_path,'thermistorstringtabelLIS122.xls'));
sensor_id=dummy(:,2);
I = find(isnan(sensor_id)==1);
sensor_id(I) = [];


calib_path=fullfile('/media/Seagate1TB/NIOZ/calibrationMatlab/');
calib_name='S_ijkbad_2013_05_08';
load(fullfile(calib_path,calib_name));
calib_path=fullfile('/media/Seagate1TB/NIOZ/calibrationMatlab/');
calib_name='S_ijkbad_2013_05_16';
load(fullfile(calib_path,calib_name));

t=thermistor();
% filename=fullfile(cruise_path,thermistor_path,'oldmod4','0106_001.SEN.dec.bin');
filename=fullfile(cruise_path,thermistor_path,'mod4c','0933_002.f7f');
t=append_data(t,filename,'nodata');
t_for_range=ref_date(t,datenum(2001,01,01));

%%%%ROC12: diurnal (tidal)!
dag1=281.142;
dag2=290.071;
ra=[6.1 9.1];
skipdata=100;
% dag1=281.8;%I(a0)...284-285!;front 282-283
% dag2=282.8;
% ra=[6.5 9];
% skipdata=20;
% dag1=282.35;%Ia...284-285!;front 282-283
% dag2=282.45;
% ra=[6.5 9];
% skipdata=3;
% dag1=282.356;%Iaa...284-285!;front 282-283
% dag2=282.371;
% ra=[7.73 8.7];
% skipdata=0;
% dag1=284;%II...284-285!;front 282-283
% dag2=285;
% ra=[6.5 9];
% skipdata=50;
% dag1=284.05;%IIa...284-285!;front 282-283
% dag2=284.25;
% ra=[7.0 9.1];
% skipdata=5;
% dag1=284.174;%IIaa...284-285!;front 282-283
% dag2=284.219;
% ra=[8.03 9.06];
% skipdata=1;
% dag1=284.35;%IIb...284-285!;front 282-283
% dag2=284.52;
% ra=[7 8.2];
% skipdata=5;
% dag1=284.425;%IIbb...284-285!;front 282-283
% dag2=284.49;
% ra=[6.92 8.0];
% skipdata=1;
% dag1=284.427;%IIbbb...284-285!;front 282-283
% dag2=284.45;
% ra=[7.35 8.05];
% skipdata=0;
% dag1=285.6;%III
% dag2=286.6;
% ra=[6.0 9];
% skipdata=20;
% dag1=286.002;%IIIaa...sink-front
% dag2=286.02;
% ra=[6.8 8.3];
% skipdata=0;
% $$$ dag1=286.5;%IIIb...front
% $$$ dag2=286.6;
% $$$ ra=[6.6 7.5];
% $$$ skipdata=3;
% dag1=286.493;%IVbb...284-285!;front 282-283
% dag2=286.55;
% ra=[6.8 7.6];
% skipdata=1;
% dag1=288.493;%Vbb...284-285!;front 282-283
% dag2=288.55;
% ra=[7.2 8];
% skipdata=1;
cdag1=281.142;
cdag2=290.071;
skipcal=1000;

% thres=0.01;
% degr=1;
% % thres=0.0017;
% % degr=4;
thres=0.0017;
degr=5;
cont=20;%...1/cont degr
range_wanted=[cdag1 cdag2];
range=range_estimator(t_for_range,range_wanted,'yeardays')
skip=skipcal;

s=string();
for p=1:length(sensor_id)
  clear t 
  strnumber=num2str(sensor_id(p),'%04d');
  if sensor_id(p)>250
      for k=1:4
          filename_tmp=fullfile(cruise_path,thermistor_path,...
                                'mod4c',[strnumber '_' num2str(k,'%03d'),'.f7f']);
          d=dir(filename_tmp);
          if ~isempty(d) && d.bytes>250e6
              filename=filename_tmp;
          end
      end
  else
      filename=fullfile(cruise_path,thermistor_path,'oldmod4',[strnumber '_001.SEN.dec.bin']);
  end
  %   display(filename);
  t=thermistor();
  if exist(filename,'file') && ~ismember(sensor_id(p),[199 198 197])

      t=append_data(t,filename,'skip',skip,'range',range,'noaccel','patch');
      %           t=append_data(t,filename,'skip',skip,'noaccel','patch');
      %           t=append_data(t,filename,'skip',skip,'noaccel');
  end
  s=append_therm(s,t); % append empty t() if not existing
end
s=ref_date(s,datenum(2001,1,1));

%%
load(fullfile(calib_path,calib_name));
s=import_calib(s,S_ijkbad_2013_05_08,S_ijkbad_2013_05_16);

a=get_value(s,'calib','start_date','id');
f=zeros(1,length(a.calib));ff=zeros(1,length(a.calib));
for k=1:length(a.calib) 
    f(k)=isempty(a.calib{k}.temperature.p);
    ff(k)=isempty(a.start_date{k}.raw.value);
end
index_nocalib=setdiff(find(f),find(ff));
a.id(index_nocalib);
%%

s=import_calib(s,index_nocalib,S_ijkbad_2013_05_08,'average');

%%
s=auto_scale(s,'plot','threshold',thres,'method','polyval','degree',degr);%bathcalib

%%
s2=auto_scale2(s,[13 15 65 85 110],'n_iter',1)

%%
s_calib=remove_data(s2);

%%
%read short period
range_wanted=[dag1 dag2];
range=range_estimator(t_for_range,range_wanted,'yeardays')
skip=skipdata;
s=string();
for p=1:length(sensor_id)
  clear t
  strnumber=num2str(sensor_id(p),'%04d');
  if sensor_id(p)>250
      for k=1:4
          filename_tmp=fullfile(cruise_path,thermistor_path,...
              'mod4c',[strnumber '_' num2str(k,'%03d'),'.f7f']);
          d=dir(filename_tmp);
          if ~isempty(d) && d.bytes>250e6
              filename=filename_tmp;
          end
      end
  else
      filename=fullfile(cruise_path,thermistor_path,'oldmod4',[strnumber '_001.SEN.dec.bin']);
  end
%   display(filename);

  t=thermistor();
      if exist(filename,'file') && ~ismember(sensor_id(p),[199 198 197])

          t=append_data(t,filename,'skip',skip,'range',range,'noaccel','patch');
%           t=append_data(t,filename,'skip',skip,'noaccel','patch');
%           t=append_data(t,filename,'skip',skip,'noaccel');
      end
  s=append_therm(s,t);
end
s=ref_date(s,datenum(2001,1,1));
s=import_calib(s,s_calib);

%%

[tijd,Torig]=string2mat(s,'calib','yearday');
[tijd,TT]=string2mat(s,'calib','yearday');
[mT,nT]=size(TT);
dz=zeros(length(sensor_id),1);
dz(1:51)=-[1:51]*0.6+mean(pA(3:nA))/1.01+120+0.6;
dz(52:length(sensor_id))=-[1:length(sensor_id)-51]*1.0+dz(51);
sa=35.3*ones(mT,nT);
dplct10=gsw_CT_from_t(sa,TT,dz*1.01-10.1325);
TT=dplct10;
nf1=1;
nf2=nT;
% midT=mean(TT(:,2:nf2-1)');
% adiab=1.86876e-4;
% adiabp=adiab*[1:103]+midT(1);
% stabtrend=(midT(mT)-adiabp(mT)-(midT(1)-adiabp(1)))/(mT-1);
% figure(1);clf;plot(midT,dz,adiabp,dz,'r')
% pause
T1=TT;
for p=[22 85]
    T1(p,:)=mean(TT(p-1:2:p+1,:));
end
for p=[13 15 48 56 58 88 126 128 130]
    T1(p,:)=mean(TT(p-1:2:p+1,:));
end
T0=T1;
jst=4;
for p=[98]%
    for jj=1:jst
        T0(p+jj-1,:)=T0(p-1,:)*(jst+1-jj)/(jst+1)+T0(p+jst,:)*jj/(jst+1);
    end
end
jst=5;
for p=[24]%
    for jj=1:jst
        T0(p+jj-1,:)=T0(p-1,:)*(jst+1-jj)/(jst+1)+T0(p+jst,:)*jj/(jst+1);
    end
end
jst=7;
for p=[90]%
    for jj=1:jst
        T0(p+jj-1,:)=T0(p-1,:)*(jst+1-jj)/(jst+1)+T0(p+jst,:)*jj/(jst+1);
    end
end

% figure;imagesc(tijd,-dz,T0);axis xy;caxis(ra);colormap('darkjet');colorbar
% break
%%

% pause

drdT=0.157*1e-2;
T=T0;
% T=flipdim(T0,1);
[mt,nt]=size(T);
nf1=1;
nf2=nt;
% % for ii=1:mt
% %     for jj=1:nt
% %        T0(ii,jj)=gsw_CT_from_t(35,T(ii,jj),dz(ii));
% %     end
% % end
% % T=T0;
deltaz=1.0;
[Tsort,I]=sort(T,1);
[tmp,J]=sort(I,1);
[time,Z]=meshgrid([1:size(T,2)],[1:size(T,1)]);
%   ThScale=(I-Z)*deltaz;
[dTdt,dTTdz]=gradient(Tsort,1,deltaz);
ddz=-gradient(dz);
for ii=1:nt
    dTdz(:,ii)=dTTdz(:,ii)./ddz;
    ThScale(:,ii)=(I(:,ii)-Z(:,ii)).*ddz;
end
kz=0.128*ThScale.^2.*(drdT*dTdz).^0.5;
flux=kz.*dTdz;
Nbv=(drdT*dTdz).^0.5;
nu=1.e-6;
eps=0.64*ThScale.^2.*Nbv.^3;
keyboard
kzhi=1.6*abs(ThScale).*(nu*Nbv).^0.5;
epsthres=100*nu*Nbv.^2;
mtijd1=tijd(nf1:nf2);
mstq1=mean(Nbv(2:mt-1,nf1:nf2).^2);
mth1=mean(kz(2:mt-1,nf1:nf2));
meps=mean(eps(2:mt-1,nf1:nf2));
mfl=mean(1e3*kz(2:mt-1,nf1:nf2).*Nbv(2:mt-1,nf1:nf2).^2);
mflkz=mfl./mean(Nbv(2:mt-1,nf1:nf2).^2)./1e3;

az=isnan(mth1);
bz=isnan(meps);
for ii=1:nt
    if az(ii)==1
        mth1(ii)=1e-10;
    end
    if bz(ii)==1
        meps(ii)=1e-10;
    end
end
shortp=2*pi./max(max(Nbv))
mmNbv=mean(mstq1).^0.5
mmkz=mean(mth1)
mmeps=mean(meps)
%  
% [R,C]=find(isnan(T));
% nf1=max(C(C<size(T,2)/2))+1;
% nf2=min(C(C>size(T,2)/2))-1;
% nf1=1;
% nf2=n;
xt1(1)=tijd(nf1)/2+tijd(nf2)/2;
xt1(2)=xt1(1)+shortp/86400;
yt1(1)=dz(mt-2);
yt1(2)=dz(mt-2);
ptr(1)=180+06/24+15/1440;
ptr(2)=180+06/24+15/1440;
ptt(1)=0.154;
ptt(2)=1.75;
xw1(1)=ptr(1)-5/1440;
xw1(2)=ptr(1)+5/1440;
yw1(1)=ptt(1)-0.03;
yw1(2)=ptt(1)-0.03;
xw2(1)=ptr(1)-5/1440;
xw2(2)=ptr(1)+5/1440;g
yw2(1)=ptt(2)-0.1;
yw2(2)=ptt(2)-0.1;
x2(1)=180.3330;
x2(2)=180.3370;
y2(1)=0.08;
y2(2)=0.08;
y3=y2;
x3(1)=180.2870;
x3(2)=180.2930;
% ra=[5.2 6.7]; 
rb=[-4 -1.9]; 
rc=[-50 50]; 
re=[-10 -4]; 
ira=[tijd(nf1) tijd(nf2) dz(length(sensor_id)) dz(1)+7];
ird=[tijd(nf1) tijd(nf2) -5 0];
dsplit=dag2-dag1;
px1=[dag1 dag1+dsplit/3 dag1+dsplit*2/3 dag2];
px0=[' ',' ',' ',' '];
% px1=[99.5 100 100.5 101 101.5 102 102.5 103 103.5 104];
% px0=[' ',' ',' ',' '];
figure(4)
clf
p1=[.15,.82,.7,.15];
axes('position',p1)
imagesc(tijd,dz,Torig,ra);colorbar;colormap('jet')
hold on
hn=plot(xt1,yt1,'m');
set(hn,'LineWidth',3)
hold off
axis(ira)
set(gca,'Fontsize',9)
ylabel('-z (m)')
set(gca,'XTick',px1);
set(gca,'XTickLabel',px0);
% lk=text(xt1(1)+0.001,yt1(1)-7,'600 s');
lk=text(xt1(1)+0.001,yt1(1)-7,num2str(shortp,'%4.0f'));
set(lk,'Fontsize',10,'Color','m')
lk=text(tijd(nf2)+0.75,915,'T (^oC)');
set(lk,'Fontsize',12)
lk=text(104.1,945,'a');
set(lk,'Fontsize',14)
% lk=text(tijda(na2)+0.00,0.1,'a');
% set(lk,'Fontsize',14)

p1=[.15,.66,.7,.15];
axes('position',p1)
imagesc(tijd,dz,Tsort,ra-0.11);colorbar;colormap('jet')
hold on
hn=plot(xt1,yt1,'m');
set(hn,'LineWidth',3)
hold off
axis(ira)
set(gca,'Fontsize',9)
ylabel('-z (m)')
set(gca,'XTick',px1);
set(gca,'XTickLabel',px0);
lk=text(tijd(nf2)+0.75,915,'\Theta (^oC)');
set(lk,'Fontsize',12)
lk=text(104.1,945,'b');
set(lk,'Fontsize',14)

p1=[.15,.5,.7,.15];
axes('position',p1)
imagesc(tijd,dz,log10(Nbv),rb);colorbar;
axis(ira)
set(gca,'Fontsize',9)
ylabel('-z (m)')
set(gca,'XTick',px1);
set(gca,'XTickLabel',px0);
lk=text(tijd(nf2)+0.75,915,'log(N) (s^-^1)');
set(lk,'Fontsize',12)
lk=text(104.1,945,'c');
set(lk,'Fontsize',14)

p1=[.15,.34,.7,.15];
axes('position',p1)
imagesc(tijd,dz,ThScale,rc);colorbar;
axis(ira)
set(gca,'Fontsize',9)
ylabel('-z (m)');xlabel('(yearday)')
set(gca,'XTick',px1);
set(gca,'XTickLabel',px0);
lk=text(tijd(nf2)+0.75,915,'d (m)');
set(lk,'Fontsize',12)
lk=text(104.1,945,'d');
set(lk,'Fontsize',14)

py1=[-4 -3 -2 -1 0];
py0=[' ',' ',' ',' '];
p1=[.15,.24,.5755,.09];
axes('position',p1)
% hn=plot(mtijd1,log10(mflkz),'r',mtijd1,log10(mth1),'k');
hn=plot(mtijd1,log10(mth1),'k',mtijd1,log10(meps)+4,'r');
set(hn,'LineWidth',0.4)
% hold on
% hn=plot(mtijd1,log10(6e-9./mstq1),'m');
% set(hn,'LineWidth',1.0)
% hold off
axis(ird)
set(gca,'Fontsize',8)
% ylabel('log(K_z)')
set(gca,'YAxisLocation','Right')
set(gca,'XTick',px1);
set(gca,'XTickLabel',px0);
set(gca,'YTick',py1);
set(gca,'YTickLabel',py1);
lk=text(tijd(nf2)+0.2,-2.5,'log(<K_z>) (m^2 s^-^1)');
set(lk,'Fontsize',11)
lk=text(103.8,-1,'e');
set(lk,'Fontsize',14)

p1=[.15,.08,.7,.15];
axes('position',p1)
hmm=colormap;
hmm(1,:)=[1 1 1];
imagesc(tijd,dz,log10(eps),re);colorbar;
colormap(hmm)
axis(ira)
set(gca,'Fontsize',9)
ylabel('-z (m)');xlabel('(yearday)')
set(gca,'XTick',px1);
set(gca,'XTickLabel',px1);
lk=text(tijd(nf2)+0.75,915,'log(\epsilon) (W kg^-^1)');
set(lk,'Fontsize',11)
lk=text(104.1,945,'f');
set(lk,'Fontsize',14)

print -djpeg I:\projects\ROCLIS12\mooring\matlab\ROC12epsdetIIIb

break
irk=[7.4 12 -dz(1) -dz(140)];
figure(1)
clf
k1=1;
Td=zeros(140,1);
Tu=zeros(140,1);
for ii=1:140
    Td(ii)=T(140-ii+1,k1+ii-1);
    Tu(ii)=T(ii,k1+ii-1);
end
subplot(1,2,1)
plot(T(:,k1),-dz,'k-',T(:,k1+140),-dz,'k-.',Td(140:-1:1),-dz,'r',Tu,-dz,'r--')
subplot(1,2,2);%buoy scale >336 s
for ii=k1:2:k1+340
    plot(T(:,ii)+0.01*(ii-k1),-dz,'k')
    axis(irk)
    hold on
end
hold off
print -djpeg ROC12buoyTest1

figure(2)
clf
ff=find(tijd<282.3585);
[kn,km]=size(ff);
k1=ff(km);
Td=zeros(140,1);
Tu=zeros(140,1);
for ii=1:140
    Td(ii)=T(140-ii+1,k1+ii-1);
    Tu(ii)=T(ii,k1+ii-1);
end
subplot(1,2,1)
plot(T(:,k1),-dz,'k-',T(:,k1+140),-dz,'k-.',Td(140:-1:1),-dz,'r',Tu,-dz,'r--')
subplot(1,2,2);%buoy scale >336 s-
for ii=k1:2:k1+340
    plot(T(:,ii)+0.01*(ii-k1),-dz,'k')
    axis(irk)
    hold on
end
hold off
print -djpeg ROC12buoyTest2

figure(3)
clf
ff=find(tijd<282.417);
[kn,km]=size(ff);
k1=ff(km);
Td=zeros(140,1);
Tu=zeros(140,1);
for ii=1:140
    Td(ii)=T(140-ii+1,k1+ii-1);
    Tu(ii)=T(ii,k1+ii-1);
end
subplot(1,2,1)
plot(T(:,k1),-dz,'k-',T(:,k1+140),-dz,'k-.',Td(140:-1:1),-dz,'r',Tu,-dz,'r--')
subplot(1,2,2);%buoy scale >336 s
for ii=k1:2:k1+340
    plot(T(:,ii)+0.01*(ii-k1),-dz,'k')
    axis(irk)
    hold on
end
hold off
print -djpeg ROC12buoyTest3


