% Basic script used to generate transect project of glider data.
%  Deployment -> M299 - "Nexos01"
close all
clear all

%% Load data file (originaly glider_contours2.m)
matFile = '~/data/GliderData/SEA003/20161118/matlab/SEA003_20161118_l2.mat'
output = glider_process_socib(matFile, 10, 325, 1);
names = fieldnames(output);
for i=1:length(names)
    eval([names{i} '=output.' names{i} ';']);
end

%% By-pass existing calibration   <------------ ***
% (ideally, netcdf files should be re-generated with new calib.)
NAP_calib = 0.3154; % N=5
NAP_blank = 0.0684;
PHE_calib = 3.0869;
PHE_blank = -0.0177;
%I = find(TRYru>.14); % <---- Try to flag bad counts See how I can do it better!
%TRYru(I) = NaN;
PHE = ( PHEru - PHE_blank)./PHE_calib;
NAP = ( TRYru - NAP_blank)./NAP_calib;


%% Some parameters (see notes in GreenLabBook - 7March2017)
WP = ...
    [ 60+45.0/60 3+12.0/60
      60+55.0/60 4+22.0/60
      60+46.8/60 3+24.6/60
      60+41.0/60 3+47.0/60
      60+52.0/60 3+35.2/60];
WP_names =  {'L1', 'L2', 'T1', 'T2' 'T3'};
timeLims = []; % empty = whole dataset


% 1. L1-L2 transect
origin = [WP(1,1) WP(1,2)]; % (Lat,Lon)
target = [WP(2,1) WP(2,2)];
timeLims = [timeVec(1) datenum(2016,11,21,6,0,0)]; % First transect                                          
%timeLims = [datenum(2016,11,21,5,59,0) datenum(2016,11,24,3,24,0)]; % Return transect


% 2. T1-T2 (TrollA-TrollB) 
% $$$ origin = [WP(3,1) WP(3,2)]; % (Lat,Lon)
% $$$ target = [WP(4,1) WP(4,2)];
% $$$ timeLims = [timeVec(225) timeVec(260)]; % 1st transect                                           
% $$$ timeLims = [timeVec(324) timeVec(349)]; % 2nd transect                                           
% $$$ timeLims = [timeVec(349) timeVec(379)]; % 3rd transect                                           
% $$$ timeLims = [timeVec(379) timeVec(407)]; % 4th transect                                           
% $$$ timeLims = [timeVec(407) timeVec(447)]; % 5th transect                                           

% For montage, see ~/research/MIO/M312-Troll01/allPortions/README



%% Compute along-transect distance
if isempty(timeLims)
    [theIndexCTD, xVecCTD] = glider_transect_projection(lonVecCTD, latVecCTD, origin, target);
    [theIndex, xVec] = glider_transect_projection(lonVec, latVec, origin, target);
else
    [theIndexCTD, xVecCTD] = glider_transect_projection(lonVecCTD, latVecCTD, origin, target, timeVecCTD, timeLims);
    [theIndex, xVec] = glider_transect_projection(lonVec, latVec, origin, target, timeVec, timeLims);
end


%% plot script

%glider_transect_projection_plottroll01
%glider_transect_projection_plottroll01_MS
