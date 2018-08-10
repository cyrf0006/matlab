function S2 = shear2hst(adcpFile, zVec, timeVec, timeSmooth)

% function velocity2hst(adcpFile, zVec, timeVec)
%
% usage ex: 
%  shear2hst('/media/Seagate1TB/NIOZ/thermistordata/ROC12/roc12.mat', zVec, timeVec,[])   

    
%% Load and get the relevant variables
load(adcpFile);

N = [SerNmmpersec/1000]'; % in m/s now
E = [SerEmmpersec/1000]';
W = [SerVmmpersec/1000]';

% build depth matrix (time variable)
zAdcp = AnDepthmm/1000; % instrument depth
zMatrix = nan(size(N));
dz = RDIBinSize;

% remove values at wrong depth (while it was lowering?)
I = abs(zAdcp - nanmean(zAdcp))>10; % flag if depth is 10m from mean
E(:,I) = NaN;
N(:,I) = NaN;
W(:,I) = NaN;

% remove 1st line
E(1,:) = NaN;
N(1,:) = NaN;
W(1,:) = NaN;

for i = 1:length(zAdcp)    
    zMatrix(:,i) = [SerBins*RDIBinSize+zAdcp(i)+RDIBin1Mid]';
end

% time vector
timeAdcp = [datenum(SerYear+2000, SerMon, SerDay, SerHour, SerMin, SerSec)]';
clear Ser* An*

%% Compute shear
if ~isempty(timeSmooth)
    disp('must average, not implemented')
    return
end

[dEx, dEz] = gradient(E);
[dNx, dNz] = gradient(N);
S2 = (dEz./dz).^2 + (dNz./dz).^2;

%% Restrict shear to desired range (to speed up following)
I = find(timeAdcp>=timeVec(1) & timeAdcp<=timeVec(end));
S2 = S2(:,I);
zMatrix = zMatrix(:,I);
timeAdcp = timeAdcp(I);
J = find(zMatrix<min(zVec) | zMatrix>max(zVec));
S2(J) = NaN;


%% Average/interp to thermistor depth
S2z = nan(length(zVec), length(timeAdcp));
for i = 1:length(timeAdcp)
    S2z(:,i) = interp1(zMatrix(:,i), S2(:,i), abs(zVec));
end

%% Average/interp to thermistor resolution
if diff(timeAdcp(1:2)) > diff(timeVec(1:2)) % adcp velocity must be interpolated
    S2 = interp1(timeAdcp, S2z', timeVec);
    S2 = S2';
else % adcp vel must be averaged
    S2 = nan(size(S2z,1), length(timeVec))
    dt = diff(timeVec(1:2));
    disp('NEVER BEEN CHECK! [K]')
    keyboard
    for i = 1:length(timeVec)
        I = find(timeAdcp>=timeVec-dt/2 & timeAdcp<timeVec+dt/2);
        S2(:,i) = nanmean(S2z(:,I), 2);
    end
end


