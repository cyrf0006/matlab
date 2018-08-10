% Basic script used to generate transect project of glider data.
%  Deployment -> M299 - "Nexos01"
clear all
close all

%% other method

matFile = '~/Data/GliderData/SEA003/20161011/matlab/SEA003_20161011_l2.mat';
load(matFile)
F255 = sauma02.data_raw.fluorescence_255_360;
F270 = sauma02.data_raw.fluorescence_270_340;
F255m = sauma02.data_raw.fluorescence_monitoring_255_360;
F270m = sauma02.data_raw.fluorescence_monitoring_270_340;
zVec = sauma02.data_raw.depth;
timeVec = sauma02.data_raw.time;
timeVec = posixtime2utc(timeVec);
chl = sauma02.data_raw.chlorophyll;
cdom = sauma02.data_raw.cdom;
bb = sauma02.data_raw.backscatter_700;
S = sauma02.data_raw.salinity;
T = sauma02.data_raw.temperature;

%I = find(zVec<=10 & ~isnan(F255));
I = find(~isnan(F255));
F255 = F255(I);
F270 = F270(I);
F255m = F255m(I);
F270m = F270m(I);
zVec = zVec(I);
timeVec1 = timeVec(I);

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
TRY_calib = sauma02.deployment.calibration.MFL.TRY_std;
PHE_calib = sauma02.deployment.calibration.MFL.PHE_spf;
NAP_calib = sauma02.deployment.calibration.MFL.NAP_spf;
TRY_blank = sauma02.deployment.calibration.MFL.TRY_blank;
PHE_blank = sauma02.deployment.calibration.MFL.PHE_blank;
NAP_blank = sauma02.deployment.calibration.MFL.NAP_blank;
DARK = sauma02.deployment.calibration.MFL.DARK;
TRY = ( ((F270-DARK)./(F270m-DARK)) - TRY_blank)./TRY_calib;
PHE = ( ((F255-DARK)./(F255m-DARK)) - PHE_blank)./PHE_calib;
NAP = ( ((F270-DARK)./(F270m-DARK)) - NAP_blank)./NAP_calib;
PHEru = (F255-DARK)./(F255m-DARK);
TRYru = ((F270-DARK)./(F270m-DARK));
I = find(timeVec1>=datenum(2016,10,11) & timeVec1<=datenum(2016,10,11));
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
xlabel('11 Oct. 2016')



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

I = find(timeVec1>=datenum(2016,10,11, 8,10,0) & timeVec1<=datenum(2016,10,11,8,53,0));
[PHEd, spike] = despike(PHE, 3, 1/240, 2, 1);
[TRYd, spike] = despike(TRY, 3, 1/240, 2, 1);

% $$$ 
% $$$ figure(11)
% $$$ clf
% $$$ plot(timeVec1(I), PHE(I)*1000, 'k')
% $$$ hold on
% $$$ plot(timeVec1(I), PHEd(I)*1000, 'r')
% $$$ plot(timeVec1(I), runmean(PHEd(I),15)*1000, 'm', 'lineWidth', 2)
% $$$ datetick
% $$$ ylabel('Phe-like (ng L^{-1})', 'fontSize', FS, 'fontWeight', 'bold')    
% $$$ xlabel('6 July 2016', 'fontSize', FS, 'fontWeight', 'bold')
% $$$ set(gca, 'ygrid', 'on')
% $$$ title('Ancien Chapeau', 'fontSize', FS, 'fontWeight', 'bold')
% $$$ ylim([-20 60])
% $$$ legend('raw', 'dspike', 'smooth', 'location', 'northWest')
% $$$ print('-dpng', '-r300', 'phe_oldCap.png')
% $$$ 
% $$$ figure(12)
% $$$ clf
% $$$ plot(timeVec1(I), TRY(I), 'k')
% $$$ hold on
% $$$ plot(timeVec1(I), TRYd(I), 'r')
% $$$ plot(timeVec1(I), runmean(TRYd(I),15), 'm', 'lineWidth', 2)
% $$$ datetick
% $$$ ylabel('Try-like (\mug L^{-1})', 'fontSize', FS, 'fontWeight', 'bold')    
% $$$ xlabel('6 July 2016', 'fontSize', FS, 'fontWeight', 'bold')
% $$$ set(gca, 'ygrid', 'on')
% $$$ title('Ancien Chapeau', 'fontSize', FS, 'fontWeight', 'bold')
% $$$ ylim([-5 35])
% $$$ legend('raw', 'dspike', 'smooth', 'location', 'northWest')
% $$$ print('-dpng', '-r300', 'try_oldCap.png')
% $$$ 
% $$$ I = find(timeVec1>=datenum(2016,07,06, 7,22,0) & timeVec1<=datenum(2016,07,06,8,16,0));

figure(13)
clf
plot(timeVec1(I), PHE(I)*1000, 'k')
hold on
plot(timeVec1(I), runmean(PHE(I),15)*1000, 'm', 'lineWidth', 2)
datetick
ylabel('Phe-like (ng L^{-1})', 'fontSize', FS, 'fontWeight', 'bold')    
xlabel('11 Oct. 2016', 'fontSize', FS, 'fontWeight', 'bold')
set(gca, 'ygrid', 'on')
xlim([timeVec1(I(1)) timeVec1(I(end))])
%title('Nouveau Chapeau', 'fontSize', FS, 'fontWeight', 'bold')
legend('raw', 'smooth', 'location', 'northWest')
print('-dpng', '-r300', 'phe_sauma02.png')

figure(14)
clf
plot(timeVec1(I), TRY(I), 'k')
hold on
plot(timeVec1(I), runmean(TRY(I),15), 'm', 'lineWidth', 2)
datetick
ylabel('Try-like (\mug L^{-1})', 'fontSize', FS, 'fontWeight', 'bold')    
xlabel('11 Oct. 2016', 'fontSize', FS, 'fontWeight', 'bold')
set(gca, 'ygrid', 'on')
xlim([timeVec1(I(1)) timeVec1(I(end))])
%title('Nouveau Chapeau', 'fontSize', FS, 'fontWeight', 'bold')
legend('raw', 'smooth', 'location', 'northWest')
print('-dpng', '-r300', 'try_sauma02.png')
