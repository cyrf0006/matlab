function Ri = hst_Ri(adcpFile, N2, zVec, timeVec, varargin)%zVec, timeVec, timeSmooth)

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

%% Compute shear (maybe to identify a length scale)
[dUx, dUz] = gradient(U);
dz = zVecAdcp(2)-zVecAdcp(1);
S2 = (dUz./dz).^2;


%% Average/interp buoyancy freq to shear resolution
if length(timeAdcp)>1
    dtAdcp = timeAdcp(2)-timeAdcp(1);
else
    dtAdcp = newTimebin;
end
if length(zVecAdcp)>1
    dzAdcp = zVecAdcp(2)-zVecAdcp(1); 
else
   zVecAdcp = newZbin;
end

N2itp = nan(length(zVecAdcp), length(timeAdcp));
for i = 1:length(timeAdcp)
    for j = 1:length(zVecAdcp)
        I = find(timeVec>=timeAdcp(i)-dtAdcp/2 & timeVec<=timeAdcp(i)+dtAdcp/2);
        J = find(zVec>=zVecAdcp(j)-dzAdcp/2 & zVec<=zVecAdcp(j)+dzAdcp/2);
        N2itp(j,i) = nanmean(nanmean(N2(J,I)));
    end
end


Ri.Ri = N2itp./S2;
Ri.timeVec = timeAdcp;
Ri.zVec = zVecAdcp;
Ri.N2 = N2itp;
Ri.S2 = S2;



