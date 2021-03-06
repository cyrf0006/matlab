% Basic script used to generate transect project of glider data.
%  
close all
clear all

%% Load data file (originaly glider_contours2.m)
matFile = '~/data/gliders_data/SEA003/20180503/matlab/SEA003_20180503_l2.mat'
output = glider_process_socib(matFile, 0, 650, 1);
names = fieldnames(output);
for i=1:length(names)
    eval([names{i} '=output.' names{i} ';']);
end

%% By-pass existing calibration   <------------ ***
% (ideally, netcdf files should be re-generated with new calib.)
% $$$ NAP_calib = 0.3154; % N=5
% $$$ NAP_blank = 0.0684;
% $$$ PHE_calib = 3.0869;
% $$$ PHE_blank = -0.0177;
% $$$ %I = find(TRYru>.14); % <---- Try to flag bad counts See how I can do it better!
% $$$ %TRYru(I) = NaN;
% $$$ PHE = ( PHEru - PHE_blank)./PHE_calib;
% $$$ NAP = ( TRYru - NAP_blank)./NAP_calib;


%% Some parameters (see notes in GreenLabBook - 7March2017)
% $$$ WP = ...
% $$$     [ 60+45.0/60 3+12.0/60
% $$$       60+55.0/60 4+22.0/60
% $$$       60+46.8/60 3+24.6/60
% $$$       60+41.0/60 3+47.0/60
% $$$       60+52.0/60 3+35.2/60];
% $$$ WP_names =  {'L1', 'L2', 'T1', 'T2' 'T3'};
% $$$ timeLims = []; % empty = whole dataset


% 1. L1-L2 transect
origin = [38.94253, 2.96987]; % (Lat,Lon)
target = [37.88908, 2.64376];

%timeLims = [datenum(2018,05,04,22,0,0) datenum(2018,05,10,0,0,0)];% First transect
timeLims = [datenum(2018,05,09,23,30,0) datenum(2018,05,13,19,0,0)]; % Return transect

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

%glider_transect_projection_plotpreswot_zoom
%glider_transect_projection_plotpreswot

