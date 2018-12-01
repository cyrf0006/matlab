% Basic script used to generate transect project of glider data.
%  
close all
clear all

%% Load data file (originaly glider_contours2.m)
matFile = '~/data/GliderData/SEA022/20181106/matlab/SEA022_20181106_l2.mat'
output = glider_process_socib(matFile, 0, 350, 1);
names = fieldnames(output);
for i=1:length(names)
    eval([names{i} '=output.' names{i} ';']);
end

% BB transect
origin = [48.733, -52.967]; % (Lat,Lon)
target = [50.32, -47.947];
timeLims = [];
timeLims = [datenum(2018,11,06,15,0,0) datenum(2018,11,18,15,0,0)]; % First transect
                                                                    %timeLims = [datenum(2018,11,18,15,0,0) datenum(2018,11,27,15,0,0)]; % Return transect 

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
