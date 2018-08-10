function Re = hst_Re(adcpFile, zVec, maxdVec, timeVec, varargin)

% function velocity2hst(adcpFile, varargin)
%
% usage ex: 
%  buoy2shear('/media/Seagate1TB/NIOZ/thermistordata/ROC12/roc12.mat',zVecReg, timeSensor1)   


%% Varargin test
if isempty(varargin)==1
    newZbin = [];
    newTimebin = [];
elseif size(varargin,2) == 2
    newZbin = strcmp(varargin{1}, 'zbin');
    newTimebin = strcmp(varargin{1}, 'timebin');
elseif size(varargin,2) == 4
    newZbin = strcmp(varargin{1}, 'zbin') + strcmp(varargin{3}, 'zbin');
    newTimebin = strcmp(varargin{1}, 'timebin') + strcmp(varargin{3}, 'timebin');
else
    disp('Wrong input... try "help hst_Re"')
    return
end

% Find the position of the 'clockMismatch' argument within the varargin cell.
% The following argument is the actual closk offset (in days).
if newZbin == 1
    I = find(strcmp(varargin, 'zbin') == 1);
    newZbin = varargin{I+1};
    myNewZbin = 1;
else
    myNewZbin = 0;
end    
if newTimebin == 1
    I = find(strcmp(varargin, 'timebin') == 1);
    newTimebin = varargin{I+1};
    myNewTimebin = 1;
else
    myNewTimebin = 0;
end    

    
%% Load and get the relevant variables
load(adcpFile);

N = [SerNmmpersec/1000]'; % in m/s now
E = [SerEmmpersec/1000]';
%W = [SerVmmpersec/1000]';
U = sqrt(N.^2 + E.^2);

% build depth matrix (time variable)
zAdcp = AnDepthmm/1000; % instrument depth
zMatrix = nan(size(N));
dz = RDIBinSize;

% remove values at wrong depth (while it was lowering?)
I = abs(zAdcp - nanmean(zAdcp))>10; % flag if depth is 10m from mean
U(:,I) = NaN;
% remove 1st line
U(1,:) = NaN;

for i = 1:length(zAdcp)    
    zMatrix(:,i) = [SerBins*RDIBinSize+zAdcp(i)+RDIBin1Mid]';
end

% time vector
timeAdcp = [datenum(SerYear+2000, SerMon, SerDay, SerHour, SerMin, SerSec)]';
clear Ser* An*


%% Restrict shear to desired range (to speed up following)
I = find(timeAdcp>=timeVec(1) & timeAdcp<=timeVec(end));
U = U(:,I);
zMatrix = zMatrix(:,I);
timeAdcp = timeAdcp(I);
J = find(zMatrix<min(zVec) | zMatrix>max(zVec));
U(J) = NaN;

zVecAdcp = nanmean(zMatrix,2);

%% If needed, vertically average ADCP data
if myNewZbin
    zVecNew = zVecAdcp(1)+newZbin/2:newZbin:zVecAdcp(end);
    newU = nan(length(zVecNew), length(timeAdcp));
    for i = 1:length(zVecNew)
        I = find(zVecAdcp>=zVecNew(i)-newZbin/2 & zVecAdcp<=zVecNew(i)+newZbin/2);
        newU(i,:) = nanmean(U(I,:),1);
    end   
    % Rename new matrices
    zVecAdcp = zVecNew;
    U = newU;
end


%% If needed, horizontally average ADCP data
if myNewTimebin
    timeAdcpNew = timeAdcp(1)+newTimebin/2:newTimebin:timeVec(end);
    newU = nan(length(zVecAdcp), length(timeAdcpNew));
    for i = 1:length(timeAdcpNew)
        I = find(timeAdcp>=timeAdcpNew(i)-newTimebin/2 & timeAdcp<=timeAdcpNew(i)+newTimebin/2);
        newU(:,i) = nanmean(U(:,I),2);
    end   
    % Rename new matrices
    timeAdcp = timeAdcpNew;
    U = newU;
end


%% Reynolds number calculation
%h = abs(zVecAdcp(end)-zVecAdcp(1));
h = max(maxdVec); % Maximum ovt in the window
nu = 1e-6;
Re = max(U).*h/nu;





