clear all
close all

%% Few parameters
hdfFile = 'myFilename.h5';
TSrelFile = 'temp_rho_relFile.mat';
TSrelFile = 'temp_rho_relFile_3rd.mat';
nSkip = 11;
smoothingWindow = 20; % minutes (for S2, N2 smoothing)
timeWindowMinutes = 60;
windowSize = timeWindowMinutes/1440; % in days
g = 9.81;
%dz = .6;
timeAveHist = [-.5:.2:.5];
totalDepth = 919;


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
%timeSensor = n+datenum(2012, 1, 1); % THIS MUST BE GENERALIZED!!!!!!!!(EX:2012)
clear n
disp('  -> done!')


%% Prepare data
% make sure matrix is flip correctly
[zVec, I] = sort(zVec);
T = T(I,:);

% reduce number of pts according to nskip
T = T(:,1:nSkip:end);
timeVec = timeSensor(1:nSkip:end);
I = find(T<-10);
T(I) = NaN;


% subsets:
% $$$ I = find(timeVec >= datenum(2012,10,8,12,25,0) & timeVec <= datenum(2012,10,8,13,25,0));%ovt
% $$$ I = find(timeVec >= datenum(2012,10,11,11,25,0) & timeVec <= datenum(2012,10,11,12,25,0));
% $$$ I = find(timeVec >= datenum(2012,10,18,19,25,0) & timeVec <= datenum(2012,10,8,20,25,0));%ovt
% $$$ I = find(timeVec >= datenum(2012,10,10,12,25,0) & timeVec <= datenum(2012,10,10,13,25,0));%ovt
% $$$ 
% $$$ I = find(timeVec >= datenum(2012,10,8,0,0,0) & timeVec <= datenum(2012,10,10,0,0,0));
I = find(timeVec >= min(timeVec) & timeVec <= max(timeVec));

T1 = T(:,I);
timeVec = timeVec(I);

%% From T to sigma
% to potential temperature
pt = gsw_pt0_from_t(T1*0+35,T1,zVec);

% to density
load(TSrelFile);
rho1 = polyval(p, T1)+1000;


%% N2 calculation
%[windowRhoSort, I] = sort(rho1,1);
%[dRx, dRz] = gradient(windowRhoSort);
fid = fopen('N2plotinfo.mat');
if fid == -1
    dz = gradient(zVec);
    N2 = nan(length(zVec), length(timeVec));
    I = find(~isnan(nanmean(rho1,2)));
    disp('  N2 calculation loop...')
    tic
    for i = 1:length(timeVec)      
        rhoVec = rho1(:,i);
        I = find(~isnan(rhoVec));
        rhoItp = interp1(zVec(I), rhoVec(I), zVec);    
        rho0 = nanmean(rhoVec);
        rhoSort = sort(rhoItp);
        dRz = gradient(rhoSort);
        dz = gradient(zVec);
        N2(:,i) = g./rho0.*dRz./dz;
        %    N2(:,i) = g./rho0.*dRz(:,i)./dz;
    end
    toc
    I = find(isnan(nanmean(T,2))); 
    N2(I,:) = NaN;  
    save N2plotinfo.mat N2
else
    load N2plotinfo.mat
    if length(N2(:)) ~= length(T(:));
        disp(['size(N2) doest fit temperature, please delete N2plotinfo.mat and rerun']);
        break
    end
end


%% read and prepare ADCP data
[U, E, N, W] = velocity2hst('./roc12.mat', zVec, timeVec) ;
[U, V] = rotate_vecd(E,N,-30);
%S2 = shear2hst('./roc12.mat', zVec, timeVec,[]);

Ri = hst_Ri('./roc12.mat',N2,zVec, timeVec);
% N2 and S2 inside Ri

%% read epsilon (saved in APEF_vs_Thorpe, manually!)
eps = load('epsilon_Thorpe_wholeTM.mat');
epsMat = eps.epsMat;
epstime = eps.timeSensor1;

%% Phase averaging
time2 = time2maxV('./roc12.mat',timeVec); 
time2_eps = time2maxV('./roc12.mat',epstime); 
time2_Ri = time2maxV('./roc12.mat',Ri.timeVec); 


% Average temperature field relative to max vel
dtAve = .05;
timeAve = [-.5:dtAve:.5];
tempAve = nan(length(zVec), length(timeAve));
UAve = nan(length(zVec), length(timeAve));
VAve = nan(length(zVec), length(timeAve));
S2Ave = nan(length(Ri.zVec), length(timeAve));
N2Ave = nan(length(Ri.zVec), length(timeAve));
RiAve = nan(length(Ri.zVec), length(timeAve));
epsAve = nan(length(eps.zVecReg), length(timeAve));
for i = 1:length(timeAve)
   I = find(time2>=timeAve(i)-dtAve & time2<=timeAve(i)+dtAve);
   tempVec = nanmean(T1(:,I),2);
   UAve(:,i) = nanmean(U(:,I),2);
   VAve(:,i) = nanmean(V(:,I),2);
   %S2Ave(:,i) = nanmean(S2(:,I),2);
   %N2Vec = nanmean(N2(:,I),2);
   %RiVec = nanmean(Ri(:,I),2);

   % interpolate missing sensors (in T)
   I = ~isnan(tempVec);
   tempVec = interp1(zVec(I), tempVec(I), zVec);   
   tempAve(:,i) = tempVec;
% $$$    N2Vec = interp1(zVec(I), N2Vec(I), zVec);   
% $$$    N2Ave(:,i) = N2Vec;
   
% $$$    RiVec = interp1(zVec(I), RiVec(I), zVec);   
% $$$    RiAve(:,i) = RiVec;

   % Ri
   I = find(time2_Ri>=timeAve(i)-dtAve & time2_Ri<=timeAve(i)+dtAve);
   RiVec = nanmean(Ri.Ri(:,I),2);
   S2Vec = nanmean(Ri.S2(:,I),2);
   N2Vec = nanmean(Ri.N2(:,I),2);
   I = isinf(RiVec);
   RiVec(I) = NaN;
   S2Vec(I) = NaN;
   N2Vec(I) = NaN;

   I = find(~isnan(RiVec));
   if sum(I) > 2
       RiVec = interp1(Ri.zVec(I), RiVec(I), Ri.zVec);   
       RiAve(:,i) = RiVec;
       S2Vec = interp1(Ri.zVec(I), S2Vec(I), Ri.zVec);   
       S2Ave(:,i) = S2Vec;
       N2Vec = interp1(Ri.zVec(I), N2Vec(I), Ri.zVec);   
       N2Ave(:,i) = N2Vec;
   end

   % epsilon
   I = find(time2_eps>=timeAve(i)-dtAve & time2_eps<=timeAve(i)+dtAve);
   epsAve(:,i) = nanmean(epsMat(:,I),2);
end

% $$$ % replace NaNs in N2
% $$$ I = find(isnan(nanmean(T1,2)));
% $$$ N2Ave(I,:) = NaN;


I = find(isnan(nanmean(UAve,2)));
I = [I; I(end)+1; I(end)+2; I(end)+3; size(UAve,1)];
UAve(I,:) = NaN;
VAve(I,:) = NaN;
I = find(sum(isnan(S2Ave),2) ~= 0);
%I = [I; I(end)+1; I(end)+2; I(end)+3; size(S2Ave,1)];
S2Ave(I,:) = NaN;
N2Ave(I,:) = NaN;

%tidal_average_plot
tidal_average_plot2
%tidal_average_plot3

keyboard



% Find phase of an event
myTime = datenum(2012,10,8,13,0,0)-.5;
myTime = datenum(2012,10,11,12,0,0)-.5;
myTime = datenum(2012,10,10,20,0,0)-.5;
[Y, I] = min(abs(timeVec - myTime));
time2(I)
