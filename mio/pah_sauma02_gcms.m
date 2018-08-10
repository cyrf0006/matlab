% Basic script used to generate transect project of glider data.
%  Deployment -> M299 - "Nexos01"
clear all
close all

%% Glider data
matFile = '~/data/GliderData/SEA003/20161011/matlab/SEA003_20161011_l2.mat';
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

% ---- Updated calibration (pah_est paper) ---- %
%NAP_calib = 0.1838;
%NAP_blank = 0.1795;
NAP_calib = 0.3154; % N=5
NAP_blank = 0.0684;
%NAP_calib = 0.3511; %N=4
%NAP_blank = 0.0531;
PHE_calib = 3.0869;
PHE_blank = -0.0177;
% --------------------------------------------- %
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
p
PHE_is = (PHEru - p(2))./p(1);
p = polyfit(NAPHS,naph_is,1);
p
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
rigs = 0.05; % very right of figure
tops = 0.05; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% $$$ % *************************************************************** %
subplot(211)
plot(timeVec1(I), NAP(I)*1000, 'color', [.5 .5 .5])
hold on
plot(timeVec1(I), runmean(NAP(I),20)*1000, 'k', 'lineWidth', 2)
plot(time_gcms, NAPHS, 'kp', 'markerSize', 10, 'lineWidth', 2)
ylim([0 850])
datetick
ylabel('Naph-like (ng L^{-1})', 'fontSize', FS, 'fontWeight', 'bold')    
set(gca, 'ygrid', 'on')
set(gca, 'xticklabel', [])
xlim([timeVec1(I(1)) timeVec1(I(end))])
text(time_gcms(1), NAPHS(1), 'S01  ', 'fontSize', FS, 'fontWeight', 'bold', 'verticalAlignment', 'bottom', 'horizontalAlignment', 'right')
text(time_gcms(2), NAPHS(2), 'S02  ', 'fontSize', FS, 'fontWeight', 'bold', 'verticalAlignment', 'bottom', 'horizontalAlignment', 'right')
text(time_gcms(3), NAPHS(3), 'S03  ', 'fontSize', FS, 'fontWeight', 'bold', 'verticalAlignment', 'bottom', 'horizontalAlignment', 'right')
text(time_gcms(4), NAPHS(4), 'S04  ', 'fontSize', FS, 'fontWeight', 'bold', 'verticalAlignment', 'bottom', 'horizontalAlignment', 'right')
text(time_gcms(5), NAPHS(5), 'S05  ', 'fontSize', FS, 'fontWeight', 'bold', 'verticalAlignment', 'bottom', 'horizontalAlignment', 'right')
%legend('MFL raw', 'MFL smooth', 'GC-MS', 'location', 'northEast')
adjust_space

subplot(212)
plot(timeVec1(I), PHE(I)*1000, 'color', [.5 .5 .5])
hold on
plot(timeVec1(I), runmean(PHE(I),20)*1000, 'k', 'lineWidth', 2)
plot(time_gcms, PHES, 'kp', 'markerSize', 10, 'lineWidth', 2)
ylim([0 70])
datetick
ylabel('Phe-like (ng L^{-1})', 'fontSize', FS, 'fontWeight', 'bold')    
xlabel('11 Oct. 2016', 'fontSize', FS, 'fontWeight', 'bold')
set(gca, 'ygrid', 'on')
xlim([timeVec1(I(1)) timeVec1(I(end))])
legend('MFL raw', 'MFL smooth', 'GC-MS', 'location', 'northEast')
adjust_space

print('-dpng', '-r300', 'sauma02_gcms.png')

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
rigs = 0.05; % very right of figure
tops = 0.05; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% $$$ % *************************************************************** %

subplot(211)
plot(timeVec1(I), NAP_is(I), 'color', [.5 .5 .5])
hold on
plot(timeVec1(I), runmean(NAP_is(I),20), 'k', 'lineWidth', 2)
plot(time_gcms, NAPHS, 'kp', 'markerSize', 10, 'lineWidth', 2)
datetick
ylabel('Naph-like (ng L^{-1})', 'fontSize', FS, 'fontWeight', 'bold')    
set(gca, 'ygrid', 'on')
set(gca, 'xticklabel', [])
xlim([timeVec1(I(1)) timeVec1(I(end))])
text(time_gcms(1), NAPHS(1), 'S01  ', 'fontSize', FS, 'fontWeight', 'bold', 'verticalAlignment', 'bottom', 'horizontalAlignment', 'right')
text(time_gcms(2), NAPHS(2), 'S02  ', 'fontSize', FS, 'fontWeight', 'bold', 'verticalAlignment', 'bottom', 'horizontalAlignment', 'right')
text(time_gcms(3), NAPHS(3), 'S03  ', 'fontSize', FS, 'fontWeight', 'bold', 'verticalAlignment', 'bottom', 'horizontalAlignment', 'right')
text(time_gcms(4), NAPHS(4), 'S04  ', 'fontSize', FS, 'fontWeight', 'bold', 'verticalAlignment', 'bottom', 'horizontalAlignment', 'right')
text(time_gcms(5), NAPHS(5), 'S05  ', 'fontSize', FS, 'fontWeight', 'bold', 'verticalAlignment', 'bottom', 'horizontalAlignment', 'right')
%legend('MFL raw', 'MFL smooth', 'GC-MS', 'location', 'northEast')
adjust_space

subplot(212)
plot(timeVec1(I), PHE_is(I), 'color', [.5 .5 .5])
hold on
plot(timeVec1(I), runmean(PHE_is(I),20), 'k', 'lineWidth', 2)
plot(time_gcms, PHES, 'kp', 'markerSize', 10, 'lineWidth', 2)
datetick
ylabel('Phe-like (ng L^{-1})', 'fontSize', FS, 'fontWeight', 'bold')    
xlabel('11 Oct. 2016', 'fontSize', FS, 'fontWeight', 'bold')
set(gca, 'ygrid', 'on')
xlim([timeVec1(I(1)) timeVec1(I(end))])
legend('MFL raw', 'MFL smooth', 'GC-MS', 'location', 'northEast')
adjust_space

print('-dpng', '-r300', 'sauma02_is_gcms.png')


figure(3)
clf
% $$$ % *********************** Adjust_space.m ************************ %
% $$$ % Fields required by the function adjust_space.m. Please fill every
% $$$ % of the following and call "adjust_space" in the script whenever
% $$$ % you want. Do not touch four last fields
ncol = 2; % no. subplot column
nrow = 2; % no. subplot row
dx = 0.08 ; % horiz. space between subplots
dy = 0.09; % vert. space between subplots
lefs = 0.08; % very left of figure
rigs = 0.05; % very right of figure
tops = 0.05; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% $$$ % *************************************************************** %

subplot(221)
plot(timeVec1(I), NAP_is(I), 'color', [.5 .5 .5])
hold on
plot(timeVec1(I), runmean(NAP_is(I),20), 'k', 'lineWidth', 2)
plot(time_gcms, NAPHS, 'kp', 'markerSize', 10, 'lineWidth', 2)
datetick
ylabel('Naph-like (ng L^{-1})', 'fontSize', FS, 'fontWeight', 'bold')    
set(gca, 'ygrid', 'on')
set(gca, 'xticklabel', [])
xlim([timeVec1(I(1)) timeVec1(I(end))])
text(time_gcms(1), NAPHS(1), 'S01  ', 'fontSize', FS, 'fontWeight', 'bold', 'verticalAlignment', 'bottom', 'horizontalAlignment', 'right')
text(time_gcms(2), NAPHS(2), 'S02  ', 'fontSize', FS, 'fontWeight', 'bold', 'verticalAlignment', 'bottom', 'horizontalAlignment', 'right')
text(time_gcms(3), NAPHS(3), 'S03  ', 'fontSize', FS, 'fontWeight', 'bold', 'verticalAlignment', 'bottom', 'horizontalAlignment', 'right')
text(time_gcms(4), NAPHS(4), 'S04  ', 'fontSize', FS, 'fontWeight', 'bold', 'verticalAlignment', 'bottom', 'horizontalAlignment', 'right')
text(time_gcms(5), NAPHS(5), 'S05  ', 'fontSize', FS, 'fontWeight', 'bold', 'verticalAlignment', 'bottom', 'horizontalAlignment', 'right')
text(timeVec1(I(1)), 86, '  a', 'fontSize', FS, 'fontWeight', 'bold')
legend('MFL raw', 'MFL smooth', 'GC-MS', 'location', 'northEast')
adjust_space

subplot(222)
plot(NAPHS/1000, naph_is, 'ok')
hold on
plot(NAP_is/1000, TRYru, 'k')
xlim([0 .1])
ylim([0 .55])
xlabel('Naphs (ug L^{-1})', 'fontSize', FS, 'fontWeight', 'bold')    
ylabel('Naph-like (RU)', 'fontSize', FS, 'fontWeight', 'bold')    
set(gca, 'ygrid', 'on')
set(gca, 'xgrid', 'on')
legend('GC-MS', 'y = 5.54x - 0.084', 'location', 'northWest')
text(.089, .5, 'b', 'fontSize', FS, 'fontWeight', 'bold')
adjust_space

subplot(223)
plot(timeVec1(I), PHE_is(I), 'color', [.5 .5 .5])
hold on
plot(timeVec1(I), runmean(PHE_is(I),20), 'k', 'lineWidth', 2)
plot(time_gcms, PHES, 'kp', 'markerSize', 10, 'lineWidth', 2)
datetick
ylabel('Phe-like (ng L^{-1})', 'fontSize', FS, 'fontWeight', 'bold')    
xlabel('11 Oct. 2016', 'fontSize', FS, 'fontWeight', 'bold')
set(gca, 'ygrid', 'on')
xlim([timeVec1(I(1)) timeVec1(I(end))])
text(timeVec1(I(1)), 18, '  c', 'fontSize', FS, 'fontWeight', 'bold')
%legend('MFL raw', 'MFL smooth', 'GC-MS', 'location', 'northEast')
adjust_space


subplot(224)
plot(PHES/1000, phe_is, 'ok')
hold on
plot(PHE_is/1000, PHEru, 'k')
xlim([0 .02])
ylim([0 .25])
xlabel('Phes (ug L^{-1})', 'fontSize', FS, 'fontWeight', 'bold')    
ylabel('Phe-like (RU)', 'fontSize', FS, 'fontWeight', 'bold')    
set(gca, 'ygrid', 'on')
set(gca, 'xgrid', 'on')
legend('GC-MS', 'y = 6.70x + 0.078', 'location', 'northWest')
text(.018, .23, 'd', 'fontSize', FS, 'fontWeight', 'bold')
adjust_space


print('-dpng', '-r300', 'sauma02_is_gcms_calib.png')
