%% Script used for in situ calibration - Troll01 campaign
% 
% Beginning of script is copy-paste from glider_contours2.m

clear all
close all


%% Troll01:
structFile = '~/data/GliderData/SEA003/20161118/matlab/SEA003_20161118_l2.mat';

zMin = 0;
zMax = 400;
downcast = 1; % =1 for downcast only; =0 for all
g = 9.81;

%% Load structure
disp(sprintf('load %s', structFile));
load(structFile);
w = whos;

found = 0;
i = 1;
while found == 0
    if strcmp(w(i).class, 'struct')
        found = 1;
        structName = w(i).name;
    else
        i = i+1;
    end
end
command = sprintf('s = %s ;',structName);
eval(command);

% Dimensions vectors
lonVec = s.data_grid.longitude;
latVec = s.data_grid.latitude;
timeVec = posixtime2utc(s.data_grid.time);
zVec = s.data_grid.depth;
distanceOrigin = nan(size(timeVec));
distanceCumul = nan(size(timeVec));
distanceCumul(1) = 0;

% $$$ for i = 1:length(timeVec)
% $$$     distanceOrigin(i) = m_lldist([lonVec(1) lonVec(i)], [latVec(1) latVec(i)]);
% $$$ end
% $$$ 
% $$$ for i = 2:length(timeVec)
% $$$     distanceCumul(i) = distanceCumul(i-1) + m_lldist([lonVec(i-1) lonVec(i)], [latVec(i-1) latVec(i)]);
% $$$ end

% lat-lon for GSW toolbox
lon = nanmean(latVec);
lat = nanmean(lonVec);


%% -------------------- Data ---------------------- %%

%% SEABIRD (.25 Hz)
disp(' -> Seabird data')
T = s.data_grid.temperature';
C = s.data_grid.conductivity'*10;
S = s.data_grid.salinity';
P = s.data_grid.pressure';

if downcast == 1 % select only downcast for CTD
    Idown = find(nanmean(s.data_grid.pitch')<0);    
    T = T(:,Idown);
    C = C(:,Idown);
    P = P(:,Idown);
    distanceCumulCTD = distanceCumul(Idown);
    latVecCTD = latVec(Idown);
    lonVecCTD = lonVec(Idown);    
else    
    distanceCumulCTD = distanceCumul;
    latVecCTD = latVec;
    lonVecCTD = lonVec;
end

% restrict depth
Iz = find(zVec>zMin & zVec<zMax);
T = T(Iz,:);
C = C(Iz,:);
S = S(Iz,:);
P = P(Iz,:);
zVec = zVec(Iz);

% Cleaning (just remove NaNs for now...)
[Cc Tc Pc] = sx_cleanCTD(timeVec, zVec, C,T,P);

% GSW stuff
SP = gsw_SP_from_C(Cc,Tc,Pc);
[SA, in_ocean] = gsw_SA_from_SP(SP,Pc,lon,lat);
CT = gsw_CT_from_t(SA,Tc,Pc);
sig0 = gsw_sigma0(SA,CT);
rho = gsw_rho(SA,CT,Pc);

% N2 calculation
[rhoSort, I] = sort(rho,1);
[dRx, dRz] = gradient(rhoSort);
[dPx, dPz] = gradient(Pc);
N2 = g./1030.*dRz./dPz;


%% WL-TRIPLET
disp(' -> Wetlab triplet')
CHL = s.data_grid.chlorophyll';
BB = s.data_grid.backscatter_700';
CDOM = s.data_grid.cdom';
CHL = CHL(Iz,:);
BB = BB(Iz,:);
CDOM = CDOM(Iz,:);

% Despike/clean?
I = find(BB<0);
BB(I) = NaN;

% Quick vertical ITP
[CHL BB CDOM] = sx_cleanCTD(timeVec, zVec, CHL,BB,CDOM);


%% MiniFLuo
disp(' -> Minifluo-UV')

% !!!!!!!!!!! by pass SN for test !!!!!!!!!!
%s.deployment.calibration.SN_MFL = 10;
if s.deployment.calibration.MFL.SN <= 8 % legacy
    TRYc = s.data_grid.fluorescence_270_340'; % counts
    PHEc = s.data_grid.fluorescence_255_360'; % counts
    PHEm = s.data_grid.fluorescence_monitoring_270_340'; % counts_monitor
    TRYm = s.data_grid.fluorescence_monitoring_255_360'; % counts_monitor
else
    TRYc = s.data_grid.fluorescence_270_340'; % counts
    PHEc = s.data_grid.fluorescence_255_360'; % counts
    TRYm = s.data_grid.fluorescence_monitoring_270_340'; % counts_monitor
    PHEm = s.data_grid.fluorescence_monitoring_255_360'; % counts_monitor
end

% restrict depths
TRYc = TRYc(Iz,:);
TRYm = TRYm(Iz,:);
PHEc = PHEc(Iz,:);
PHEm = PHEm(Iz,:);

% clean (must be modified)
[TRYc, TRYm, YY] = sx_cleanCTD(timeVec, zVec, TRYc,TRYm,TRYc);
[PHEc, PHEm, YY] = sx_cleanCTD(timeVec, zVec, PHEc,PHEm,TRYc);

% get + apply calibration
TRY_calib = s.deployment.calibration.MFL.TRY_std;
PHE_calib = s.deployment.calibration.MFL.PHE_spf;
NAP_calib = s.deployment.calibration.MFL.NAP_spf;
TRY_blank = s.deployment.calibration.MFL.TRY_blank;
PHE_blank = s.deployment.calibration.MFL.PHE_blank;
NAP_blank = s.deployment.calibration.MFL.NAP_blank;
DARK = s.deployment.calibration.MFL.DARK;

% ---- Updated calibration (pah_est paper) ---- %
%NAP_calib = 0.1838;
%NAP_blank = 0.1795;
NAP_calib = 0.3154;
NAP_blank = 0.0684
PHE_calib = 3.0869;
PHE_blank = -0.0177;
% --------------------------------------------- %


% Lab Calibration
TRY_ru = ( ((TRYc-DARK)./(TRYm-DARK)) );
PHE_ru = ( ((PHEc-DARK)./(PHEm-DARK)) );
if isempty(TRY_blank) % Assume blank = 0;
    TRY = TRY_ru./TRY_calib;
    PHE = PHE_ru./PHE_calib;
    NAP = TRY_ru./NAP_calib;
else
    TRY = (TRY_ru - TRY_blank)./TRY_calib;
    PHE = (PHE_ru - PHE_blank)./PHE_calib;
    NAP = (TRY_ru - NAP_blank)./NAP_calib;
end

% Manual blank
NAP2 = (TRY_ru - .066)./NAP_calib;
PHE2 = (PHE_ru - .04)./PHE_calib;


% ---------------------------------------------------- %


% Initial profile
nCasts = 1;
NAPinit = nanmean(NAP(:,1:nCasts),2); %ng/L
PHEinit = nanmean(PHE(:,1:nCasts),2);
NAPinit2 = nanmean(NAP2(:,1:nCasts),2); %ng/L
PHEinit2 = nanmean(PHE2(:,1:nCasts),2);
NAPinit_ru = nanmean(TRY_ru(:,1:nCasts),2); %ng/L
PHEinit_ru = nanmean(PHE_ru(:,1:nCasts),2);


% In situ data
dz = 2;
Cphes = [12.6 1.9 24.0 1.0];
Cnaphs = [152.4 186.6 137 81.1];
zPAHs = [0 10 25 40];
zPAHs = [1 8 22 35];
zPAHs = [2 7 21 35];

Phe_ru_gcsm = nan(size(zPAHs));
Nap_ru_gcsm = nan(size(zPAHs));
for i = 1:length(zPAHs)
    I = find(zVec>=zPAHs(i)-2 & zVec<=zPAHs(i)+2);
    Phe_ru_gcms(i) = nanmean(PHEinit_ru(I));
    Nap_ru_gcms(i) = nanmean(NAPinit_ru(I));
end


%% plots 
FS = 14;
FS2 = 12;
figure(1)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 10 15])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 1; % no. subplot row
dx = 0.0 ; % horiz. space between subplots
dy = 0.0; % vert. space between subplots
lefs = 0.15; % very left of figure
rigs = 0.05; % very right of figure
tops = 0.02; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** 

plot(NAPinit.*1000, zVec, 'k', 'lineWidth', 2)
hold on
plot(NAPinit2.*1000, zVec, '--k', 'lineWidth', 2)
plot(Cnaphs, zPAHs, 'mp', 'markerSize', 10, 'markerFaceColor', 'm')
set(gca, 'ydir', 'reverse')
ylim([0 120])
%title('Average of first 5 glider profiles vs water sample', 'fontSize', FS, 'fontWeight', 'bold')
xlabel('C_{Naphs} (ngL^{-1})', 'fontSize', FS, 'fontWeight', 'bold')
ylabel('Depth (m)', 'fontSize', FS, 'fontWeight', 'bold')
legend('MFL fluorescence', 'Blank adjusted', 'GC/MS', 'location', 'southWest')
set(gca, 'fontSize', FS2)
adjust_space
print(gcf, '-dpng', 'Naphs_GCMS.png')

figure(2)
clf
set(gcf,'PaperUnits','centimeters','PaperPosition',[1 1 10 15])
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 1; % no. subplot row
dx = 0.0 ; % horiz. space between subplots
dy = 0.0; % vert. space between subplots
lefs = 0.15; % very left of figure
rigs = 0.05; % very right of figure
tops = 0.02; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** 
plot(PHEinit.*1000, zVec, 'k', 'lineWidth', 2)
hold on
plot(PHEinit2.*1000, zVec, '--k', 'lineWidth', 2)
plot(Cphes, zPAHs, 'mp', 'markerSize', 10, 'markerFaceColor', 'm')
set(gca, 'ydir', 'reverse')
ylim([0 120])
%title('Average of first 5 glider profiles vs water sample', 'fontSize', FS, 'fontWeight', 'bold')
xlabel('C_{Phes} (ngL^{-1})', 'fontSize', FS, 'fontWeight', 'bold')
ylabel('Depth (m)', 'fontSize', FS, 'fontWeight', 'bold')
legend('MFL fluorescence', 'Blank adjusted', 'GC/MS', 'location', 'southWest')
set(gca, 'fontSize', FS2)
adjust_space
print(gcf, '-dpng', 'Phes_GCMS.png')


%% Return MiniFluo values for comparision with GC-MS
dz=5
mfl_phe = nan(size(zPAHs))
mfl_nap = nan(size(zPAHs))
for i = 1:length(zPAHs)
    I = find(zVec>=zPAHs(i)-dz/2 & zVec<=zPAHs(i)+dz/2);
    mfl_phe(i) = nanmean(PHEinit(I));
    mfl_nap(i) = nanmean(NAPinit(I));
end
   
disp('Mfl returned concentration for Naph:')
mfl_nap*1000
disp('Mfl returned concentration for Phe:')
mfl_phe*1000