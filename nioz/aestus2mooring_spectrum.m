
figure(3)
clf
set(gcf, 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 16 18])

%% some info
zSpec = 850; %m


%% ADCP
adcpFile = '../../roc12.mat';
load(adcpFile);

N = [SerNmmpersec/1000]'; % in m/s now
E = [SerEmmpersec/1000]';
W = [SerVmmpersec/1000]';

[U, V] = rotate_vecd(E,N,-30); %<--- maybe 30deg. is not the best angle with model...

% build depth matrix (time variable)
zAdcp = AnDepthmm/1000; % instrument depth
zMatrix = nan(size(N));

% remove values at wrong depth (while it was lowering?)
I = abs(zAdcp - nanmean(zAdcp))>10; % flag if depth is 10m from mean
U(:,I) = NaN;
V(:,I) = NaN;

for i = 1:length(zAdcp)    
    zMatrix(:,i) = [SerBins*RDIBinSize+zAdcp(i)+RDIBin1Mid]';
end

% time vector
timeAdcp = [datenum(SerYear+2000, SerMon, SerDay, SerHour, SerMin, SerSec)]';
clear Ser* An*
freqAdcp = 1/((timeAdcp(2)-timeAdcp(1))*86400); 

% timeseries for spectral analysis
zAdcp = nanmean(zMatrix,2);
[Y, I] = min(abs(zAdcp-zSpec));
vAdcp = V(I,:);
uAdcp = U(I,:);
I = find(isnan(vAdcp));
vAdcp(I) = 0;
I = find(isnan(uAdcp));
uAdcp(I) = 0;

%% model
freqModel = 1/((timeVec(2)-timeVec(1))*86400); 
[Y, I] = min(abs(zuVec-zSpec));
vModel = uMat(I,:); % <--- v in mooring reference
uModel = vMat(I,:);
I = find(isnan(vModel));
vModel(I) = 0;
I = find(isnan(uModel));
uModel(I) = 0;


%% Spectral analysis
nx = max(size(vAdcp)); % Hanning
na = 1;
w = hanning(floor(nx/na));
vAdcp = detrend(vAdcp);
uAdcp = detrend(uAdcp);
[ps_vAdcp, f_vAdcp] = pwelch(vAdcp, w, 0, [], freqAdcp*86400); 
[ps_uAdcp, f_uAdcp] = pwelch(uAdcp, w, 0, [], freqAdcp*86400); 

nx = max(size(vModel)); % Hanning
na = 1;
w = hanning(floor(nx/na));
vModel = detrend(vModel);
uModel = detrend(uModel);
[ps_vModel, f_vModel] = pwelch(vModel, w, 0, [], freqModel*86400); 
[ps_uModel, f_uModel] = pwelch(uModel, w, 0, [], freqModel*86400); 

subplot(211)
loglog(f_vAdcp, ps_vAdcp, 'k', 'linewidth', 2)  
hold on
loglog(f_vModel, ps_vModel, 'r', 'linewidth', 2)  
xlim([1e-1 3e1])
ylim([1e-11 1e0])
plot_harmonics6
text(1.1e-1, 1e-1, 'cross-isobaths', 'horizontalAlignment', 'left', ...
     'fontWeight', 'bold')
set(gca, 'xticklabel', [])

subplot(212)
loglog(f_uAdcp, ps_uAdcp, 'k', 'linewidth', 2)  
hold on
loglog(f_uModel, ps_uModel, 'r', 'linewidth', 2)  
xlim([1e-1 3e1])
ylim([1e-11 1e0])
plot_harmonics6
text(1.1e-1, 1e-1, 'along-isobaths', 'horizontalAlignment', 'left', ...
     'fontWeight', 'bold')
legend('Mooring', 'Model', 'location', 'southWest')


print('-dpng', '-r300', 'spectra_MooringModel.png')
