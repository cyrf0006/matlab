function [y, spike] = sx_despikeBinCtd (data, zbin, smoothwindow, stdThresh)




% low-pass timeserie
cpm = 1/zbin;
nx = max(size(myT)); % Hanning
%timeWindow = 30*60; %sec. window
timeWindow = 10*60; %sec. window (for iow_BBL_contours.m)
na = length(myTime)./(timeWindow*fs);
w = hanning(floor(nx/na));

