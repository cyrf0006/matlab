function output = glider_process_socib(structFile, zMin, zMax, downcast)

% Function similar to glider_contours2, but now using .mat file
% generated by the SOCIB glider_toolbox as an input 
%
% In a further release, the function should be rename
% sx_"something" and use varargin as arguments.
%
% usage ex:
% >> output = glider_process_socib('~/data/gliders_data/SEA007/20151028/matlab/SEA007_20151028_l2.mat', 10, 250, 1);
% >> output = glider_process_socib('~/data/gliders_data/SEA003/20160324/matlab/SEA003_20160324_l2.mat' 10, 250, 1);
% >> output = glider_process_socib('~/data/gliders_data/SEA003/20160429/matlab/SEA003_20160429_l2.mat' 10, 250, 1);
% >> output = glider_process_socib('~/data/gliders_data/SEA003/20160429/matlab/SEA003_20160429_l2.mat' 10, 325, 1);

%% Add this in varargin eventually:
% $$$ zMin = 10;
% $$$ zMax = 250;
% $$$ downcast = 1; % =1 for downcast only; =0 for all
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
for i = 1:length(timeVec)
    distanceOrigin(i) = m_lldist([lonVec(1) lonVec(i)], [latVec(1) latVec(i)]);
end
for i = 2:length(timeVec)
    distanceCumul(i) = distanceCumul(i-1) + m_lldist([lonVec(i-1) lonVec(i)], [latVec(i-1) latVec(i)]);
end

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
O = s.data_grid.oxygen_frequency';
if downcast == 1 % select only downcast for CTD
    Idown = find(nanmean(s.data_grid.pitch')<0);    
    T = T(:,Idown);
    C = C(:,Idown);
    P = P(:,Idown);
    O = O(:,Idown);
    distanceCumulCTD = distanceCumul(Idown);
    latVecCTD = latVec(Idown);
    lonVecCTD = lonVec(Idown);    
    timeVecCTD = timeVec(Idown);
else    
    distanceCumulCTD = distanceCumul;
    latVecCTD = latVec;
    lonVecCTD = lonVec;
    timeVecCTD = timeVec;
end

% restrict depth
Iz = find(zVec>zMin & zVec<zMax);
T = T(Iz,:);
C = C(Iz,:);
S = S(Iz,:);
P = P(Iz,:);
O = O(Iz,:);
zVec = zVec(Iz);

% Cleaning (just remove NaNs for now...)
[Cc Tc Pc] = sx_cleanCTD(timeVecCTD, zVec, C,T,P);

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


% Oxygen concentration
[Oc, ~, ~] = sx_cleanCTD(timeVecCTD, zVec, O,T,P);
Soc = s.deployment.calibration.SBE43F.soc;
Foffset = s.deployment.calibration.SBE43F.foffset;
A = s.deployment.calibration.SBE43F.a;
B = s.deployment.calibration.SBE43F.b;
C = s.deployment.calibration.SBE43F.c;
E = s.deployment.calibration.SBE43F.e;

K = Tc + 273.15; % Temp in Kelvin
O2sol = gsw_O2sol(SA,CT,Pc,lon,lat); %umol/Kg
O2 =  Soc.*(Oc+Foffset).*(1.0+A.*CT + B.*CT.^2 + C.*CT.^3).*O2sol.*exp(E.*Pc./K); 


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
% Check version of MFL
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
    
TRY_ru = ( ((TRYc-DARK)./(TRYm-DARK)) );
PHE_ru = ( ((PHEc-DARK)./(PHEm-DARK)) );

if isempty(TRY_blank) % Assume blank = 0;
    TRY = ( ((TRYc-DARK)./(TRYm-DARK)) )./TRY_calib;
    PHE = ( ((PHEc-DARK)./(PHEm-DARK)) )./PHE_calib;
    NAP = ( ((TRYc-DARK)./(TRYm-DARK)) )./NAP_calib;
else
    TRY = ( ((TRYc-DARK)./(TRYm-DARK)) - TRY_blank)./TRY_calib;
    PHE = ( ((PHEc-DARK)./(PHEm-DARK)) - PHE_blank)./PHE_calib;
    NAP = ( ((TRYc-DARK)./(TRYm-DARK)) - NAP_blank)./NAP_calib;
end


output = struct();
output.timeVec = timeVec;
output.timeVecCTD = timeVecCTD;
output.zVec = zVec;
output.zMin = zMin;
output.zMax = zMax;
output.latVec = latVec;
output.latVecCTD = latVecCTD;
output.lonVec = lonVec;
output.lonVecCTD = lonVecCTD;
output.CT = CT;
output.SA = SA;
output.sig0 = sig0;
output.N2 = N2;
output.O2 = O2;
output.CHL = CHL;
output.BB = BB;
output.CDOM = CDOM;
output.TRY = TRY;
output.PHE = PHE;
output.NAP = NAP;
output.TRYru = TRY_ru;
output.PHEru = PHE_ru;