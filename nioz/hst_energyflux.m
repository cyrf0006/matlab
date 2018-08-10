% Once ready, make a function
clear all
close all

%% MAYBE WE SHOULD ERASE THIS FUNCTION.... 6 May 2015

%% Few parameters
hdfFile = 'myFilename.h5';
TSrelFile = 'temp_rho_relFile_3rd.mat';
nSkip = 0;
smoothingWindow = 20; % minutes (for S2, N2 smoothing)
timeWindowMinutes = 60;
windowSize = timeWindowMinutes/1440; % in days
g = 9.81;
dz = .6;
totalDepth = 919;
V = 5:.05:10;
myRi = []; myGa = [];


%% read whole HDF file
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
disp('  -> done!')

%% Prepare data
% make sure matrix is flip correctly
[zVec, I] = sort(zVec);
T = T(I,:);

% reduce number of pts according to nskip
if nSkip ~=0
    T = T(:,1:nSkip:end);
    timeSensor = timeSensor(1:nSkip:end);
end
I = find(T<-10);
T(I) = NaN;


% subsets:
%I1 = find(timeSensor >= datenum(2012,10,8,12,25,0)-.5 & timeSensor<= datenum(2012,10,8,13,25,0)-.5);%ovt
%I1 = find(timeSensor >= datenum(2012,10,11,11,25,0)-.5 & timeSensor <= datenum(2012,10,11,12,25,0)-.5);%cvt
%I1 = find(timeSensor >= datenum(2012,10,8,19,25,0)-.5 & timeSensor <= datenum(2012,10,8,20,25,0)-.5);%low1
%I1 = find(timeSensor >= datenum(2012,10,10,12,25,0)-.5 & timeSensor <= datenum(2012,10,10,13,25,0)-.5);%low2
I1 = find(timeSensor >= datenum(2012,10,10,20,25,0) & timeSensor <= datenum(2012,10,10,24,25,0));%Bore 4h
% $$$ I1 = find(timeSensor >= datenum(2012,10,10,20,25,0) & timeSensor <= datenum(2012,10,10,26,25,0));%Bore2
% $$$ I1 = find(timeSensor >= datenum(2012,10,10,20,25,0) & timeSensor <= datenum(2012,10,11,9,25,0));%Bore2
%I1 = find(timeSensor >= datenum(2012,10,10, 5,0,0) & timeSensor <= datenum(2012,10,10,8,0,0));%double-breaking
%I1 = find(timeSensor >= datenum(2012,10,10, 6,20,0) & timeSensor <= datenum(2012,10,10,7,20,0));%double-breaking(zoom)

%I1 = find(timeSensor >= datenum(2012,10,10,23,0,0) & timeSensor <= datenum(2012,10,11,0,0,0));%cvt

%I1 = find(timeSensor >= datenum(2012,10,10,12,25,0)-.5+ii/24 & timeSensor ...
%          <= datenum(2012,10,10,13,25,0)-.5+ii/24);%low2

%I1 = find(timeSensor >= datenum(2012,10,8,0,0,0)-.5 & timeSensor <= datenum(2012,10,10,0,0,0)-.5);%48h
%I1 = find(timeSensor >= min(timeSensor) & timeSensor <= max(timeSensor));

T1 = T(:,I1);
timeSensor1 = timeSensor(I1);


%% From T to sigma
% to potential temperature
pt = gsw_pt0_from_t(T1*0+35,T1,zVec);

% to density
load(TSrelFile);
rho1 = polyval(p, T1)+1000;


%% Vertical interpolation for missing sensors and regular grid
disp('Vectical inpterpolation for missing sensors...')
zVecReg = zVec(1):dz:zVec(end);
Ritp1 = [];
Titp = [];
for i = 1:length(timeSensor1)
    rVec = rho1(:,i);
    II = find(~isnan(rVec) & rVec~=0);
    Ritp1 = [Ritp1 [interp1(zVec(II), rVec(II), zVecReg)]'];
    
    TVec = T1(:,i);
    II = find(~isnan(TVec) & TVec~=0);
    Titp = [Titp [interp1(zVec(II), TVec(II), zVecReg)]'];
end
disp('  -> done!')

%% By-pass the preceeding and load existing file
load apef_thorpe_wholeTM.mat

%% Read ADCP data
vel =  velocity2hst('roc12.mat', zVecReg, timeSensor1);

%% Reorder vertical profiles
disp('Reordering profiles')
windowRho = Ritp1;
windowTime = timeSensor1;
H = size(windowRho,1);





keyboard




% N2 calculation
[windowRhoSort, I] = sort(windowRho,1);
[dRx, dRz] = gradient(windowRhoSort);
N2bg = g./windowRhoSort.*dRz./dz;

%% Compute Ri
%[Re, Ri] = turbulent_regime('/media/Seagate1TB/NIOZ/thermistordata/ROC12/roc12.mat',N2bg);
%Ri = hst_Ri('./roc12.mat',N2bg,zVecReg, timeSensor1, 'zbin', 10, 'timebin', .5/24);
Ri = hst_Ri('./roc12.mat',N2bg,zVecReg, timeSensor1);

APEFVec = []; 
N2Vec = []; 
S2Vec = []; 
JbVec = []; JbVec1 = []; JbVec2 = []; JbVec_eddy = [];
epsVec = []; epsVec1 = []; epsVec2 = [];
epsMat = nan(length(zVecReg), length(windowTime));
maxdVec = [];
for j = 1:length(windowTime)
    rhoVec_raw =  windowRho(:,j);
    [rhoVec_sort, I] = sort(rhoVec_raw);
    d = zVecReg - zVecReg(I);
    EP1 = sum(rhoVec_raw.*(max(zVecReg)-zVecReg')); 
    EP2 = sum(rhoVec_sort.*(max(zVecReg)-zVecReg'));
    rho0 = nanmean(windowRho(:,j));
    maxdVec = [maxdVec max(d)];
    
    % local stratification
    N2local = g./rho0.*gradient(rhoVec_sort)./gradient(zVecReg');
    N2Vec = [N2Vec nanmean(N2local)];   
    APEFVec = [APEFVec g./H/rho0.*(EP1-EP2)];
    JbVec = [JbVec g./H/rho0.*(EP1-EP2).*nanmean(sqrt(N2local))];
    JbVec_eddy = [JbVec_eddy g./H/rho0.*(EP1-EP2)./(rms(d).^(2/3).*(.64*rms(d).^2.*nanmean(sqrt(N2local)).^3).^(-1/3))];

    epsVec = [epsVec .64*rms(d).^2.*nanmean(sqrt(N2local)).^3];

    %local shear
    S2Vec = [S2Vec nanmean(S2bg(:,j))];
    
    % background stratification
    I = find(windowTime>(windowTime(j)-smoothingWindow/1440/2) & windowTime<(windowTime(j)+smoothingWindow/1440/2));
    N2 = N2bg(:,I);
    N2 = N2(:);
    N = sqrt(nanmean(N2));
    S2 = S2bg(:,I);
    S2 = S2(:);
    S = sqrt(nanmean(S2));
    JbVec1 = [JbVec1 g./H/rho0.*(EP1-EP2).*N];
    epsVec1 = [epsVec1 .64*rms(d).^2.*N.^3]; 
    
    JbVec2 = [JbVec2 g./H/rho0.*(EP1-EP2).*S];
    epsVec2 = [epsVec2 .64*rms(d).^2.*N.^3]; 
    
    epsMat(:,j) = .64*d.^2.*N.^3;
end
disp('  -> done!')
 

%% Compute Re
%Re = hst_Re('./roc12.mat',zVecReg, timeSensor1, 'zbin', 10, 'timebin', .5/24);
Re = hst_Re('./roc12.mat',zVecReg, maxdVec, timeSensor1);
Re = nanmean(Re);

%% Few plots.

APEF_vs_Thorpe_plot3
%APEF_vs_Thorpe_plot4


%end
keyboard



%% Mixing efficiency distribution
GamVec = JbVec./epsVec;
I = find(GamVec>1);
GamVec(I) = 1;
dclass = .05;
classes = dclass/2:dclass:1-dclass/2;
GamClass = nan(size(classes));
for i = 1:length(classes)
    I = find(GamVec>classes(i)-dclass/2 & GamVec<=classes(i)+dclass/2);
    GamClass(i) = length(I)./length(GamVec);
end 

figure(2)
clf
set(gcf, 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 15 12])
bar(classes, GamClass, 'k')
xlabel('\Gamma')
ylabel('f')
hold on
med = median(GamVec);
plot([med med], [0 max(GamClass)], '--k')
text(med+dclass/2, max(GamClass)-dclass/5, sprintf('med = %3.2f',med))
text(med+dclass/2, max(GamClass)-dclass/2.5, sprintf('n = %d',length(GamVec)))
xlim([0 1])
ylim([0 ceil(max(GamClass)*100)/100])
set(gca, 'fontSize', 10)
set(gcf, 'renderer', 'painters')
print('-dpng', 'GammaHistogram.png')
print('-depsc2','GammaHistogram.eps')

