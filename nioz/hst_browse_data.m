function hst_browse_data(configFile)

% usage ex: hst_browse_data('ROCconfig')

% VERIF means that the portion of script has been verified that it
% was not inducing error...
    
% read config file    
eval(configFile);

% load calibration files
for i = 1:size(calibFiles,1)
    fname = calibFiles(i,:);
    I = find(fname==' '); % remove empty space if there is  
                          %fname(I) = [];
    cal{i} = load(fname);
    load(fname)
    clear fname
end

% USe an existing list or generate it (see configFile)
if isempty(thermistor_list)
    thermList = hst_build_sensor_list(sensor_id, fileMinSize, '/media/Seagate1TB/NIOZ/thermistordata/ROC12/mod4c', '/media/Seagate1TB/NIOZ/thermistordata/ROC12/oldmod4');
    disp('*** A list have been generated ***') 
    disp('(look in ./sensors.list for overview)')
    disp('  -> and/or press key to continue')
    pause
else
    fid = fopen(thermistor_list);
    C = textscan(fid, '%s', 'delimiter', '\n'); 
end

% initialize thermistor
t=thermistor();

% This might be a OCTAVE PROB.
if size(thermList, 1) == 141
    thermList  = thermList(2:end);
end
 
t=append_data(t,char(thermList(1)),'nodata');

% Loop for calibration;
t_for_range=ref_date(t,REF_TIME);
range_wanted=[cdag1 cdag2];
range=range_estimator(t_for_range,range_wanted,'yeardays');
skip=skipcal;
s=string();
for p=1:length(sensor_id)
    clear t
    t=thermistor();
    filename = char(thermList(p));
    display(filename);
    if strcmp(filename, 'missing') ~= 1 & isempty(filename) ~= 1  &  strcmp(filename, 'missing_too_short') ~= 1
        t=append_data(t,filename,'skip',skip,'range',range,'noaccel','patch');    
    end
    s=append_therm(s,t);
end
s=ref_date(s,datenum(2001,1,1));


% Calibration   (VERIF.)
% this is weak, but we suppose whether 1,2 or 3 calibfiles
if size(cal,2) == 1
    cal1 = struct2cell(cal{1});
    s=import_calib(s, cal1{1});
elseif size(cal,2) == 2
    cal1 = struct2cell(cal{1});
    cal2 = struct2cell(cal{2});
    s=import_calib(s, cal1{1}, cal2{1});
elseif size(cal,2) == 3
    cal1 = struct2cell(cal{1});
    cal2 = struct2cell(cal{2});
    cal3 = struct2cell(cal{3});
    s=import_calib(s, cal1{1}, cal2{1}, cal3{1});
else
   disp('More than 3 calib files have been found, we use 3')
   cal1 = struct2cell(cal{1});
   cal2 = struct2cell(cal{2});
   cal3 = struct2cell(cal{3});
   s=import_calib(s, cal1{1}, cal2{1}, cal3{1});
end
% replace this line:
% $$$ s=import_calib(s,S_ijkbad_2013_05_08,S_ijkbad_2013_05_16);


% Some tricks on calibration that I don't get 
a=get_value(s,'calib','start_date','id');
f=zeros(1,length(a.calib));ff=zeros(1,length(a.calib));
for k=1:length(a.calib) 
    f(k)=isempty(a.calib{k}.temperature.p);
    ff(k)=isempty(a.start_date{k}.raw.value);
end
index_nocalib=setdiff(find(f),find(ff));
a.id(index_nocalib);
s=import_calib(s,index_nocalib,cal1{1},'average'); %(VERIF)
%s=import_calib(s,index_nocalib,S_ijkbad_2013_05_08,'average');

% HERE!!!
%[s, shift]=hst_autoscale(s,'plot','treshold',thres,'method', 'polyval','degree',degr, 'shift'); %Test
                                                                                                

% auto-scale
%s=auto_scale(s,'plot','treshold',thres,'method','polyval','degree',degr);%bathcalib
[s, shift]=auto_scale(s,'plot','treshold',thres,'method', 'polyval','degree',degr, 'shift');%bathcalib

% Autoscale 2 to correct few sensors
I = find(abs(shift)>.02); % this could be optimized
%s2=auto_scale2(s,I,'n_iter',1);                          
s2=auto_scale2(s,[13 15 65 85 110],'n_iter',1);
s3=auto_scale2(s2,[13 15],'n_iter',1); % test, not working
                                       % really... - FC


% Here you can by-pass auto_scale2:  % -FC
s_calib=remove_data(s2);
%s_calib=remove_data(s);


% Loop to build the string (shorter period):
s=string();
range_wanted=[dag1 dag2]; %read short period
range=range_estimator(t_for_range,range_wanted,'yeardays')
skip=skipdata;
s=string();
for p=1:length(sensor_id)
    clear t
    t=thermistor();
    filename = char(thermList(p));
    display(filename);
    if strcmp(filename, 'missing') ~= 1 & isempty(filename) ~= 1  &  strcmp(filename, 'missing_too_short') ~= 1
        t=append_data(t,filename,'skip',skip,'range',range,'noaccel','patch');    
    end
    s=append_therm(s,t);
end
s=ref_date(s,datenum(2001,1,1));
s=import_calib(s,s_calib);

% Build the matrix
[tijd,Torig]=string2mat(s,'calib','yearday');
TT = Torig; %[tijd,TT]=string2mat(s,'calib','yearday');
[mT,nT]=size(TT);
timeVec = tijd + datenum(YYYY, 1, 1, 0, 0, 0);

% Test coment
dz=zeros(length(sensor_id),1);
dz(1:51)=-[1:51]*0.6+mean(pA(3:nA))/1.01+120+0.6;
dz(52:length(sensor_id))=-[1:length(sensor_id)-51]*1.0+dz(51);
% $$$ sa=35.3*ones(mT,nT);
% $$$ dplct10=gsw_CT_from_t(sa,TT,dz*1.01-10.1325);
% $$$ TT=dplct10;
% $$$ nf1=1;
% $$$ nf2=nT;

% midT=mean(TT(:,2:nf2-1)');
% adiab=1.86876e-4;
% adiabp=adiab*[1:103]+midT(1);
% stabtrend=(midT(mT)-adiabp(mT)-(midT(1)-adiabp(1)))/(mT-1);
% figure(1);clf;plot(midT,dz,adiabp,dz,'r')
% pause

% Interpolation between missing sensors (weak, but works here)
tempVec = nanmean(Torig,2);
if isempty(Iremove)
    % smooth fit
    I = find(~isnan(tempVec));
    p = polyfit(I, tempVec(I), 3);
    x = [1:length(tempVec)]';
    tempFit = p(1)*x.^3 + p(2)*x.^2 + p(3)*x + p(4);
    
    % remove outliers
    sd = nanstd(abs(tempVec - tempFit));
    I = find(abs(tempVec-tempFit) > 3*sd);
    tempVec(I) = NaN;

    % Show missing values
    I = find(isnan(tempVec));
    disp([' !! No pts to remve was suggested. We suggest removing ' ...
          'these sensors: !!'])
    sprintf( '%d ', I)
else
    tempVec(Iremove) = NaN; % use existing remove list
end    

% the intperolation
I = find(~isnan(tempVec)==1);
for i = 1:size(TT, 2)
    TT(:,i) = interp1(I, TT(I,i), 1:size(TT,1));
end

if interactive
    figure(1)
    clf
    imagesc(TT)
    bool = 1;
    while bool
        R = input('Does Figure 1 looks okay? and do you want to continue (y/n)', 's');
        I = find(R==' '); % remove white space
        R(I) = [];
        if strcmp(R, 'n') == 1
            disp( '-> ABORT NOW!')
            return         
        elseif strcmp(R, 'y') == 1
            bool = 0;
        else
            disp(' please answer y or n')
        end
    end
end
   


keyboard

% figure;imagesc(tijd,-dz,T0);axis xy;caxis(ra);colormap('darkjet');colorbar
% break
%%

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
xw2(2)=ptr(1)+5/1440;
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

disp('[END]')
keyboard

%break
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
subplot(1,2,2);%buoy scale >336 s
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


