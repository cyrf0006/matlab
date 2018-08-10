function reorder2D(hdfFile, TSrelFile)

% function reorder2D(hdfFile)
%
% usage ex:
%  reorder2D('myFilename.h5', 'temp_rho_relFile.mat')
%

%% Few parameters
g = 9.81;
win = 10;
skip = 10;
% time range for processing
t1 = datenum(2012, 10, 9, 11, 0, 0);
t2 = datenum(2012, 10, 10, 12, 0, 0);
tWindow = 30*60; % averaging window in sec.
windowSize = tWindow/86400; % in days
dz = .6;


%% read HDF file
T = hdf5read(hdfFile, '/data/temperature');
T = T';
n = hdf5read(hdfFile, '/dims/time'); % raw sensor time
zVec = abs(hdf5read(hdfFile, '/dims/depth'));
t0_str = hdf5read(hdfFile, '/', 'StartDate');
samplingInterval = hdf5read(hdfFile, '/', 'SamplingInterval');
t0 = datenum(t0_str.Data, 'dd-mm-yyyy HH:MM:SS');
%timeSensor = n+datenum(2012, 1, 1); % THIS MUST BE GENERALIZED!!!!!!!!(EX:2012)
timeSensor = [t0:samplingInterval/86400:t0+samplingInterval/86400*(length(n)-1)]';
clear n

% make sure matrix is flip correctly
[zVec, I] = sort(zVec);
T = T(I,:);


%% Vertical interpolation for missing sensors and regular grid
I = find(timeSensor>t1 & timeSensor<t2);
timeVec = timeSensor(I(1:skip:end));
zVecReg = zVec(1):dz:zVec(end);
Titp = [];
for i = 1:skip:length(I)
    %    disp(sprintf('%d / %d', i, length(I)))
    TVec = T(:,I(i));
    II = find(~isnan(TVec) & TVec~=0);
    Titp = [Titp [interp1(zVec(II), TVec(II), zVecReg)]'];
end

figure(1)
clf
imagesc(timeVec, zVecReg, Titp)
set(gca, 'ydir', 'reverse')
title('Whole timeseries')
datetick
colorbar
ylabel('Depth (m)')
hold on

%% loop on 2D-windows
timeVecWindow = timeVec(1)+windowSize/2:windowSize/2:timeVec(end)-windowSize/2;
Jb = nan(length(timeVecWindow),1);
Jb2 = nan(length(timeVecWindow),1);
for i = 1:length(timeVecWindow)    
    disp(datestr(timeVecWindow(i)));
    I = find(timeVec>=timeVecWindow(i)-windowSize/2 & timeVec<timeVecWindow(i)+windowSize/2);

    if mod(i,2)==0
        plot([timeVecWindow(i)+windowSize/2 timeVecWindow(i)+windowSize/2], [min(zVec) max(zVec)], '--k')
        plot([timeVecWindow(i)-windowSize/2 timeVecWindow(i)-windowSize/2], [min(zVec) max(zVec)], '--k')
    else
        plot([timeVecWindow(i)+windowSize/2 timeVecWindow(i)+windowSize/2], [min(zVec) max(zVec)], '-k')
        plot([timeVecWindow(i)-windowSize/2 timeVecWindow(i)-windowSize/2], [min(zVec) max(zVec)], '-k')    
    end
    
    windowTime = timeVec(I);
    windowT = Titp(:,I);
    windowRho = hst_temp2sal(windowT, TSrelFile); 

    [xmatrix, zMatrix] = meshgrid(windowTime, zVecReg);
    rhoList = windowRho(:);
    [Y, I] = sort(rhoList);
    windowSort = reshape(Y, size(windowRho'));
    windowSort = windowSort';
    
% $$$     % uncomment to plot
% $$$     figure(2)
% $$$     clf
% $$$     imagesc(windowRho)
% $$$     
% $$$     figure(3)
% $$$     clf
% $$$     imagesc(windowSort)
% $$$     keyboard

    % remove steps in matrix:
    windowSort2 = nan(size(windowSort));
    for j = 1:size(windowSort,1)
        windowSort2(j,:) = nanmean(windowSort(j,:));
    end

    [dx,dz] = gradient(zMatrix);
    H = size(windowSort,1);
    L = size(windowSort,2);
    rho0 = nanmean(rhoList);

    EP1 = sum(sum(windowRho.*(max(zVecReg)-zMatrix)));
    EP2 = sum(sum(windowSort.*(max(zVecReg)-zMatrix)));
    EP3 = sum(sum(windowSort2.*(max(zVecReg)-zMatrix)));

    APEF = g./H./L./rho0.*(EP1-EP3);
    N2 = g/rho0.*gradient(windowSort(:,1))./gradient(zVecReg');
    Jb(i) = APEF.*nanmean(sqrt(N2));

    
    % traditional method
    APEFVec = [];
    for j = 1:length(windowTime)
        rhoVec_raw =  windowRho(:,j);
        rhoVec_sort = sort(rhoVec_raw);
        EP1 = sum(rhoVec_raw.*(max(zVecReg)-zVecReg')); 
        EP2 = sum(rhoVec_sort.*(max(zVecReg)-zVecReg'));
        rho0 = nanmean(windowRho(:,j));
        APEFVec = [APEFVec g./H/rho0.*(EP1-EP2)];   
    end    
    Jb2(i) = nanmean(APEFVec.*nanmean(sqrt(N2)));
    
end


figure(2)
semilogy(timeVecWindow, Jb, 'b')
hold on
semilogy(timeVecWindow, Jb2, '--b')
datetick
xlabel(datestr(timeVec(1),1))   
ylabel('J_b (m^2 s^{-3})')


keyboard


