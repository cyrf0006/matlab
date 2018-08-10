% Basic script used to generate transect project of glider data.
%  Deployment -> M299 - "Nexos01"
clear all
close all

%% Glider data
matFile = '~/Data/GliderData/SEA003/20161011/matlab/SEA003_20160706_l2.mat';
load(matFile)
F255 = sauma01.data_raw.fluorescence_255_360;
F270 = sauma01.data_raw.fluorescence_270_340;
F255m = sauma01.data_raw.fluorescence_monitoring_255_360;
F270m = sauma01.data_raw.fluorescence_monitoring_270_340;
zVec = sauma01.data_raw.depth;
timeVec = sauma01.data_raw.time;
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
I = find(timeVec1>=datenum(2016,10,11) & timeVec1<=datenum(2016,10,11));


%% GCMS
NAPHS = [37.3, 42.0, 54.6, 51.3, 76.2];
PHES = [0.7, 0.8, 2.1, 6.6, 15.5];
time_gcms = [datenum(2016,10,11,8,52,0),... 
             datenum(2016,10,11,8,42,0), ...
             datenum(2016,10,11,8,36,0), ...
             datenum(2016,10,11,8,28,0), ...
             datenum(2016,10,11,8,24,0)];

%% in situ calib
phe_is = nan(size(PHES));
naph_is = nan(size(PHES));
phe_compa = nan(size(PHES));
naph_compa = nan(size(PHES));
for i =1:length(time_gcms)
    I = find(timeVec1>=time_gcms(i)-1/1440 & timeVec1<=time_gcms(i)+1/1440);
    phe_is(i) = nanmean(PHEru(I));
    naph_is(i) = nanmean(TRYru(I));
    phe_compa(i) = nanmean(PHE(I));
    naph_compa(i) = nanmean(NAP(I));   
end
naph_compa*1000./NAPHS
phe_compa*1000./PHES


%plot(phe_is, PHES, 'k.')
p = polyfit(PHES,phe_is,1);
PHE_is = (PHEru - p(2))./p(1);
p = polyfit(NAPHS,naph_is,1);
NAP_is = (TRYru - p(2))./p(1);


%% plot

FS = 12;
I = find(timeVec1>=datenum(2016,10,11, 8,10,0) & timeVec1<=datenum(2016,10,11,8,53,0));

figure(1)
clf
% $$$ % *********************** Adjust_space.m ************************ %
% $$$ % Fields required by the function adjust_space.m. Please fill every
% $$$ % of the following and call "adjust_space" in the script whenever
% $$$ % you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 2; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.1; % very right of figure
tops = 0.05; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% $$$ % *************************************************************** %
subplot(211)
plot(timeVec1(I), PHE(I)*1000, 'k')
hold on
plot(timeVec1(I), runmean(PHE(I),15)*1000, 'm', 'lineWidth', 2)
plot(time_gcms, PHES, 'rp')
datetick
ylabel('Phe-like (ng L^{-1})', 'fontSize', FS, 'fontWeight', 'bold')    
set(gca, 'xticklabel', [])

set(gca, 'ygrid', 'on')
xlim([timeVec1(I(1)) timeVec1(I(end))])
legend('MFL raw', 'MFL smooth', 'GC-MS', 'location', 'northWest')
adjust_space

subplot(212)
plot(timeVec1(I), NAP(I)*1000, 'k')
hold on
plot(timeVec1(I), runmean(NAP(I),15)*1000, 'm', 'lineWidth', 2)
plot(time_gcms, NAPHS, 'rp')
datetick
ylabel('Naph-like (ng L^{-1})', 'fontSize', FS, 'fontWeight', 'bold')    
xlabel('11 Oct. 2016', 'fontSize', FS, 'fontWeight', 'bold')
set(gca, 'ygrid', 'on')
xlim([timeVec1(I(1)) timeVec1(I(end))])
adjust_space

print('-dpng', '-r300', 'sauma01_gcms.png')

figure(2)
clf
% $$$ % *********************** Adjust_space.m ************************ %
% $$$ % Fields required by the function adjust_space.m. Please fill every
% $$$ % of the following and call "adjust_space" in the script whenever
% $$$ % you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 2; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.1; % very right of figure
tops = 0.05; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% $$$ % *************************************************************** %

subplot(211)
plot(timeVec1(I), PHE_is(I), 'k')
hold on
plot(timeVec1(I), runmean(PHE_is(I),15), 'm', 'lineWidth', 2)
plot(time_gcms, PHES, 'rp')
datetick
ylabel('Phe-like (ng L^{-1})', 'fontSize', FS, 'fontWeight', 'bold')    
set(gca, 'xticklabel', [])

set(gca, 'ygrid', 'on')
xlim([timeVec1(I(1)) timeVec1(I(end))])
legend('MFL raw', 'MFL smooth', 'GC-MS', 'location', 'northWest')
adjust_space

subplot(212)
plot(timeVec1(I), NAP_is(I), 'k')
hold on
plot(timeVec1(I), runmean(NAP_is(I),15), 'm', 'lineWidth', 2)
plot(time_gcms, NAPHS, 'rp')
datetick
ylabel('Naph-like (ng L^{-1})', 'fontSize', FS, 'fontWeight', 'bold')    
xlabel('11 Oct. 2016', 'fontSize', FS, 'fontWeight', 'bold')
set(gca, 'ygrid', 'on')
xlim([timeVec1(I(1)) timeVec1(I(end))])
adjust_space

print('-dpng', '-r300', 'sauma01_is_gcms.png')
