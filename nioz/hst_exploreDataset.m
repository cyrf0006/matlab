function hst_exploreDataset(hdfFile, nSkip, timeWindowMinutes)

% function hst_exploreDataset(hdfFile, nSkip, timeWindowMinutes)
%
% usage ex:
%  hst_exploreDataset('myFilename.h5', 100, 120)
% 

%% Few parameters
windowSize = timeWindowMinutes/1440; % in days
dz = .6;



%% read HDF file
disp('loading data...')
T = hdf5read(hdfFile, '/data/temperature');
T = T';
n = hdf5read(hdfFile, '/dims/time'); % raw sensor time
zVec = abs(hdf5read(hdfFile, '/dims/depth'));
t0_str = hdf5read(hdfFile, '/', 'StartDate');
samplingInterval = hdf5read(hdfFile, '/', 'SamplingInterval');
t0 = datenum(t0_str.Data, 'dd-mm-yyyy HH:MM:SS');
timeSensor = [t0:samplingInterval/86400:t0+samplingInterval/86400*(length(n)-1)]';
clear n
disp('  ->done!')

%% Prepare data
% make sure matrix is flip correctly
[zVec, I] = sort(zVec);
T = T(I,:);

% reduce number of pts according to nskip
T = T(:,1:nSkip:end);
timeSensor = timeSensor(1:nSkip:end);
I = find(T==0);
T(I) = NaN;


%% Vertical interpolation for missing sensors and regular grid
% $$$ disp('Vectical inpterpolation for missing sensors...')
% $$$ zVecReg = zVec(1):dz:zVec(end);
% $$$ Titp = [];
% $$$ for i = 1:length(timeSensor)
% $$$     if mod(i, 10000) == 0        
% $$$         disp(sprintf('%d / %d', i, length(I)))
% $$$     end
% $$$     TVec = T(:,i);
% $$$     II = find(~isnan(TVec) & TVec~=0);
% $$$     Titp = [Titp [interp1(zVec(II), TVec(II), zVecReg)]'];
% $$$ end
% $$$ disp('  ->done!')

keyboard

%% Loop data and plot
% time vector for plotting
timeVec = timeSensor(1)+windowSize/2:windowSize:timeSensor(end)-windowSize/2;

for i = 1:length(timeVec)
    figure(1)
    clf
    I = find(timeSensor>=timeVec(i)-windowSize/2 &  timeSensor<=timeVec(i)+windowSize/2);
    imagesc(timeSensor(I), zVec, T(:,I));
    datetick
    xlim([timeSensor(I(1)) timeSensor(I(end))])
    colorbar
    xlabel(sprintf(datestr(timeSensor(I(1)), 1)))
    pause
end

