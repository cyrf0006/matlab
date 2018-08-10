% Basic script used to generate transect project of glider data.
%  Deployment -> M299 - "Nexos01"
clear all
%close all

%% other method

matFile = '~/Data/GliderData/SEA003/20160706/matlab/SEA003_20160706_l2.mat';
load(matFile)
F255 = sauma01.data_raw.fluorescence_255_360;
F270 = sauma01.data_raw.fluorescence_270_340;
F255m = sauma01.data_raw.fluorescence_monitoring_255_360;
F270m = sauma01.data_raw.fluorescence_monitoring_270_340;
zVec = sauma01.data_raw.depth;
timeVec = sauma01.data_raw.time;
latVec = sauma01.data_raw.latitude;
lonVec = sauma01.data_raw.longitude;
timeVec = posixtime2utc(timeVec);
chl = sauma01.data_raw.chlorophyll;
cdom = sauma01.data_raw.cdom;
bb = sauma01.data_raw.backscatter_700;
S = sauma01.data_raw.salinity;
T = sauma01.data_raw.temperature;

%I = find(zVec<=10 & ~isnan(F255));
I = find(~isnan(F255));
F255 = F255(I);
F270 = F270(I);
F255m = F255m(I);
F270m = F270m(I);
zVec = zVec(I);
timeVec1 = timeVec(I);
lonVec1 = lonVec(I);
latVec1 = latVec(I);

I = find(~isnan(chl));
bb = bb(I);
chl = chl(I);
cdom = cdom(I);
timeVec2 = timeVec(I);

I = find(~isnan(S));
S = S(I);
T = T(I);
timeVec3 = timeVec(I);



% calib
TRY_calib = sauma01.deployment.calibration.MFL.TRY_std;
PHE_calib = sauma01.deployment.calibration.MFL.PHE_spf;
NAP_calib = sauma01.deployment.calibration.MFL.NAP_spf;
TRY_blank = sauma01.deployment.calibration.MFL.TRY_blank;
PHE_blank = sauma01.deployment.calibration.MFL.PHE_blank;
NAP_blank = sauma01.deployment.calibration.MFL.NAP_blank;
DARK = sauma01.deployment.calibration.MFL.DARK;
TRY = ( ((F270-DARK)./(F270m-DARK)) - TRY_blank)./TRY_calib;
PHE = ( ((F255-DARK)./(F255m-DARK)) - PHE_blank)./PHE_calib;
NAP = ( ((F270-DARK)./(F270m-DARK)) - NAP_blank)./NAP_calib;
PHEru = (F255-DARK)./(F255m-DARK);
TRYru = ((F270-DARK)./(F270m-DARK));
I = find(timeVec1>=datenum(2016,07,06) & timeVec1<=datenum(2016,07,07));
% $$$ 
% $$$ figure(1)
% $$$ clf
% $$$ plot(timeVec1(I), PHEru(I))
% $$$ datetick
% $$$ ylabel('Phe-like (RU)')    
% $$$ xlabel('6 July 2016')
figure(1)
clf
plot(timeVec1(I), PHE(I)*1000)
datetick
ylabel('Phe-like')    
xlabel('6 July 2016')



% $$$ 
% $$$ I = find(timeVec2>=datenum(2016,07,06) & timeVec2<=datenum(2016,07,07));
% $$$ 
% $$$ figure(3)
% $$$ clf
% $$$ plot(timeVec2(I), chl(I))
% $$$ datetick
% $$$ ylabel('chl')    
% $$$ xlabel('6 July 2016')
% $$$ 
% $$$ figure(4)
% $$$ clf
% $$$ plot(timeVec2(I), log10(cdom(I)))
% $$$ datetick
% $$$ ylabel('cdom')    
% $$$ xlabel('6 July 2016')
% $$$ 
% $$$ figure(5)
% $$$ clf
% $$$ plot(timeVec2(I), log10(bb(I)))
% $$$ datetick
% $$$ ylabel('turbidity')    
% $$$ xlabel('6 July 2016')
% $$$ 
% $$$ 
% $$$ I = find(timeVec3>=datenum(2016,07,06) & timeVec3<=datenum(2016,07,07));
% $$$ figure(6)
% $$$ clf
% $$$ plot(timeVec3(I), S(I))
% $$$ datetick
% $$$ ylabel('S')    
% $$$ xlabel('6 July 2016')
% $$$ 
% $$$ figure(7)
% $$$ clf
% $$$ plot(timeVec3(I), T(I))
% $$$ datetick
% $$$ ylabel('T')    
% $$$ xlabel('6 July 2016')


%% Zoom on shorter timeseries


FS = 12;

I = find(timeVec1>=datenum(2016,07,06, 8,35,0) & timeVec1<=datenum(2016,07,06,9,14,0));
[PHEd, spike] = despike(PHE, 3, 1/240, 2, 1);
[TRYd, spike] = despike(TRY, 3, 1/240, 2, 1);


figure(11)
clf
plot(timeVec1(I), PHE(I)*1000, 'k')
hold on
plot(timeVec1(I), PHEd(I)*1000, 'r')
plot(timeVec1(I), runmean(PHEd(I),15)*1000, 'm', 'lineWidth', 2)
datetick
ylabel('Phe-like (ng L^{-1})', 'fontSize', FS, 'fontWeight', 'bold')    
xlabel('6 July 2016', 'fontSize', FS, 'fontWeight', 'bold')
set(gca, 'ygrid', 'on')
title('Ancien Chapeau', 'fontSize', FS, 'fontWeight', 'bold')
ylim([-20 60])
legend('raw', 'dspike', 'smooth', 'location', 'northWest')
print('-dpng', '-r300', 'phe_oldCap.png')

figure(12)
clf
plot(timeVec1(I), TRY(I), 'k')
hold on
plot(timeVec1(I), TRYd(I), 'r')
plot(timeVec1(I), runmean(TRYd(I),15), 'm', 'lineWidth', 2)
datetick
ylabel('Try-like (\mug L^{-1})', 'fontSize', FS, 'fontWeight', 'bold')    
xlabel('6 July 2016', 'fontSize', FS, 'fontWeight', 'bold')
set(gca, 'ygrid', 'on')
title('Ancien Chapeau', 'fontSize', FS, 'fontWeight', 'bold')
ylim([-5 35])
legend('raw', 'dspike', 'smooth', 'location', 'northWest')
print('-dpng', '-r300', 'try_oldCap.png')

I = find(timeVec1>=datenum(2016,07,06, 7,22,0) & timeVec1<=datenum(2016,07,06,8,16,0));

figure(13)
clf
plot(timeVec1(I), PHE(I)*1000, 'k')
hold on
plot(timeVec1(I), runmean(PHE(I),15)*1000, 'm', 'lineWidth', 2)
datetick
ylabel('Phe-like (ng L^{-1})', 'fontSize', FS, 'fontWeight', 'bold')    
xlabel('6 July 2016', 'fontSize', FS, 'fontWeight', 'bold')
set(gca, 'ygrid', 'on')
title('Nouveau Chapeau', 'fontSize', FS, 'fontWeight', 'bold')
legend('raw', 'smooth', 'location', 'northWest')
print('-dpng', '-r300', 'phe_newCap.png')

figure(14)
clf
plot(timeVec1(I), TRY(I), 'k')
hold on
plot(timeVec1(I), runmean(TRY(I),15), 'm', 'lineWidth', 2)
datetick
ylabel('Try-like (\mug L^{-1})', 'fontSize', FS, 'fontWeight', 'bold')    
xlabel('6 July 2016', 'fontSize', FS, 'fontWeight', 'bold')
set(gca, 'ygrid', 'on')
title('Nouveau Chapeau', 'fontSize', FS, 'fontWeight', 'bold')
legend('raw', 'smooth', 'location', 'northWest')
print('-dpng', '-r300', 'try_newCap.png')



%% In situ calibration
Inew = find(timeVec1>=datenum(2016,07,06, 7,22,0) & timeVec1<=datenum(2016,07,06,8,9,0));
Cphe = [3.2 43.3 11.4 0.7]; %not ordered
Ctime = nan(size(Cphe));
Cpheru = nan(size(Cphe));
Ctime(1) = datenum(2016,07,06, 8,24,0);
Ctime(2) = datenum(2016,07,06, 8,46,0);
Ctime(3) = datenum(2016,07,06, 8,54,0);
Ctime(4) = datenum(2016,07,06, 9,14,0);
Clat = [43.3560 43.3551 43.3555 43.3475];
Clon = [5.3155 5.3260 5.3207 5.3087];
Cdt = 60/86400; % 60 sec.


%% !!! The problem here is that I need to use hand GPS (no hit on GPS
%% since underwater!!!!)
% $$$ map_saumaty
% $$$ for i = 1:length(Inew)
% $$$     m_plot(lonVec1(Inew(i)), latVec1(Inew(i)), '.', 'color', 'r', 'markerSize', 5);
% $$$ end
% $$$ 
% $$$ for i = 1:length(Cphe)
% $$$     [Y, I] = min(abs((Clat(i) - latVec1(Inew))) + abs((Clon(i) - lonVec1(Inew))));
% $$$     datestr(timeVec1(I))
% $$$     I = find(timeVec1>=timeVec1(Inew(I))-Cdt/2 & timeVec1<=timeVec1(Inew(I))+Cdt/2);
% $$$     %PHEru(I)
% $$$     Cpheru(i) = nanmean(PHEru(I));
% $$$     Y
% $$$ end
% $$$     

% approx time
Ctime(4) = datenum(2016,07,06, 7,33,0);
Ctime(1) = datenum(2016,07,06, 7,48,0);
Ctime(3) = datenum(2016,07,06, 7,55,0);
Ctime(2) = datenum(2016,07,06, 8,5,0);
for i = 1:length(Cphe)
    I = find(timeVec1>=Ctime(i)-Cdt/2 & timeVec1<=Ctime(i)+Cdt/2);
    Cpheru(i) = nanmean(PHEru(I));
    Cnapru(i) = nanmean(TRYru(I));
end
    
p = polyfit(Cpheru, Cphe/1000, 1);
manualBlank = .49
PHE_insitu = (PHEru-manualBlank).*p(1);
xVec = .485:.001:.53;
yfit = polyval(p, xVec);

figure(9)
clf
plot(Cpheru, Cphe, '.k')                                              
hold on
plot(xVec, yfit*1000, 'r')    
ylabel('Measured concentration')
xlabel('RU counts')


I = find(timeVec1>=datenum(2016,07,06, 7,22,0) & timeVec1<=datenum(2016,07,06,8,16,0));
figure(10)
clf
plot(timeVec1(I), PHE_insitu(I)*1000, 'k')
hold on
plot(timeVec1(I), runmean(PHE_insitu(I),15)*1000, 'm', 'lineWidth', 2)
plot(Ctime, Cphe, 'rp', 'markerFaceColor', 'r', 'markerSize', 8)
datetick
ylabel('Phe-like (ng L^{-1})', 'fontSize', FS, 'fontWeight', 'bold')    
xlabel('6 July 2016', 'fontSize', FS, 'fontWeight', 'bold')
set(gca, 'ygrid', 'on')
title('in situ calibration', 'fontSize', FS, 'fontWeight', 'bold')
legend('raw', 'smooth', 'waterSamples', 'location', 'northWest')
print('-dpng', '-r300', 'phe_insitu.png')



figure(9)
clf
plot(timeVec1(I), PHE(I)*1000-192, 'k')
hold on
plot(timeVec1(I), runmean(PHE(I),15)*1000-192, 'm', 'lineWidth', 2)
plot(Ctime, Cphe, 'rp', 'markerFaceColor', 'r', 'markerSize', 8)
datetick
ylabel('Phe-like (ng L^{-1})', 'fontSize', FS, 'fontWeight', 'bold')    
xlabel('6 July 2016', 'fontSize', FS, 'fontWeight', 'bold')
set(gca, 'ygrid', 'on')
title('Nouveau Chapeau', 'fontSize', FS, 'fontWeight', 'bold')
legend('raw', 'smooth', 'location', 'northWest')
print('-dpng', '-r300', 'SiTaRaison.png')
