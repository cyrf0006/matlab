clear

%% Few parameters
hdfFile = 'myFilename.h5';
TSrelFile = 'temp_rho_relFile.mat';
nSkip = 10;
timeWindowMinutes = 60;
windowSize = timeWindowMinutes/1440; % in days
g = 9.81;

dz = .6;



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

% subsets:
I1 = find(timeSensor >= datenum(2012,10,8,12,25,0) & timeSensor <= datenum(2012,10,8,13,25,0));
I2 = find(timeSensor >= datenum(2012,10,11,10,25,0) & timeSensor <= datenum(2012,10,10,12,25,0));
%I1 = find(timeSensor >= datenum(2012,10,9,0,0,0) & timeSensor <= datenum(2012,10,15,0,0,0));

% $$$ myTimes = [1:23];
% $$$ for ii = 1:length(myTimes)
% $$$ I1 = find(timeSensor >= datenum(2012,10,9,myTimes(ii),25,0) & timeSensor <= datenum(2012,10,9,myTimes(ii)+1,25,0));


T1 = T(:,I1);
T2 = T(:,I2);
rho1 = hst_temp2sal(T1, TSrelFile)+1000; 
rho2 = hst_temp2sal(T2, TSrelFile)+1000; 
timeSensor1 = timeSensor(I1);
timeSensor2 = timeSensor(I2);



%% Vertical interpolation for missing sensors and regular grid
disp('Vectical inpterpolation for missing sensors...')
zVecReg = zVec(1):dz:zVec(end);
Ritp1 = [];
for i = 1:length(timeSensor1)
    rVec = rho1(:,i);
    II = find(~isnan(rVec) & rVec~=0);
    Ritp1 = [Ritp1 [interp1(zVec(II), rVec(II), zVecReg)]'];
end
Ritp2 = [];
for i = 1:length(timeSensor2)
    rVec = rho2(:,i);
    II = find(~isnan(rVec) & rVec~=0);
    Ritp2 = [Ritp2 [interp1(zVec(II), rVec(II), zVecReg)]'];
end
disp('  ->done!')



%% Fake-field (Rankine ovt)
load rankine_pi
fakeRho = flipud(rankine.rho);
fakeZ = min(zVec):(max(zVec)-min(zVec))./(size(fakeRho,1)-1):max(zVec);
faket = min(timeSensor1):(max(timeSensor1)-min(timeSensor1))./(size(fakeRho,2)-1):max(timeSensor1);



%% Reorder in 2D windows
% ---- Window 1 ---- %
windowRho = Ritp1;
windowTime = timeSensor1;
[xmatrix, zMatrix] = meshgrid(windowTime, zVecReg);
rhoList = windowRho(:);
zList = zMatrix(:);
[Y, I] = sort(rhoList);
d = zList - zList(I);
%d1Matrix = reshape(d1, size(windowRho'));
windowSort = reshape(Y, size(windowRho'));
windowSort = windowSort';

% remove steps in matrix:
windowSort2 = nan(size(windowSort));
for j = 1:size(windowSort,1)
    windowSort2(j,:) = nanmean(windowSort(j,:));
end
windowSort = windowSort2;

[dx,dz] = gradient(zMatrix);
H = size(windowSort,1);
L = size(windowSort,2);
rho0 = nanmean(rhoList);

EP1 = sum(sum(windowRho.*(max(zVecReg)-zMatrix)));
EP2 = sum(sum(windowSort.*(max(zVecReg)-zMatrix)));

APEF = g./H./L./rho0.*(EP1-EP2);
N2 = g/rho0.*gradient(windowSort(:,1))./gradient(zVecReg');
Jb1_2D = APEF.*nanmean(sqrt(N2));
eps1_2D = .64*rms(d).^2.*nanmean(sqrt(N2)).^3;

% traditional method
APEFVec = [];
epsVec1 = [];
for j = 1:length(windowTime)
    rhoVec_raw =  windowRho(:,j);
    [rhoVec_sort, I] = sort(rhoVec_raw);
    d = zVecReg - zVecReg(I);
    EP1 = sum(rhoVec_raw.*(max(zVecReg)-zVecReg')); 
    EP2 = sum(rhoVec_sort.*(max(zVecReg)-zVecReg'));
    rho0 = nanmean(windowRho(:,j));
    APEFVec = [APEFVec g./H/rho0.*(EP1-EP2)];
    epsVec1 = [epsVec1 .64*rms(d).*nanmean(sqrt(N2)).^3];
end    
Jb1_Vec = APEFVec.*nanmean(sqrt(N2));
Jb1_trad = nanmean(APEFVec.*nanmean(sqrt(N2)));
esp1_1D = nanmean(epsVec1);

% ---- Window 2 ---- %
windowRho = fakeRho;
windowTime = faket;
[xmatrix, zMatrix] = meshgrid(windowTime, fakeZ);
zList = zMatrix(:);
rhoList = windowRho(:);
[Y, I] = sort(rhoList);
d = zList - zList(I);
windowSort = reshape(Y, size(windowRho'));
windowSort = windowSort';

% remove steps in matrix:
windowSort2 = nan(size(windowSort));
for j = 1:size(windowSort,1)
    windowSort2(j,:) = nanmean(windowSort(j,:));
end
windowSort = windowSort2;

[dx,dz] = gradient(zMatrix);
H = size(windowSort,1);
L = size(windowSort,2);
rho0 = nanmean(rhoList);


EP1 = sum(sum(windowRho.*(max(fakeZ)-zMatrix)));
EP2 = sum(sum(windowSort.*(max(fakeZ)-zMatrix)));

APEF = g./H./L./rho0.*(EP1-EP2);
N2 = g/rho0.*gradient(windowSort(:,1))./gradient(fakeZ');
Jb2_2D = APEF.*nanmean(sqrt(N2));
eps2_2D = .64*rms(d).^2.*nanmean(sqrt(N2)).^3;

% traditional method
APEFVec2 = [];
epsVec2 = [];
for j = 1:length(windowTime)
    rhoVec_raw =  windowRho(:,j);
    [rhoVec_sort, I] = sort(rhoVec_raw);
    d = fakeZ - fakeZ(I);
    EP1 = sum(rhoVec_raw.*(max(fakeZ)-fakeZ')); 
    EP2 = sum(rhoVec_sort.*(max(fakeZ)-fakeZ'));
    rho0 = nanmean(windowRho(:,j));
    APEFVec2 = [APEFVec2 g./H/rho0.*(EP1-EP2)];   
    epsVec2 = [epsVec2 .64*rms(d).^2.*nanmean(sqrt(N2)).^3];
end    
Jb2_Vec = APEFVec2.*nanmean(sqrt(N2));
Jb2_trad = nanmean(APEFVec2.*nanmean(sqrt(N2)));
esp2_1D = nanmean(epsVec2);



%% Few plots
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 4; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.14; % very right of figure
tops = 0.03; % top of figure
bots = 0.07; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %
figure(1)
clf
set(gcf, 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 16 18])
s1 = subplot(4,1,[1 2]);
imagesc(timeSensor1, zVecReg, Ritp1-1000)
c = colorbar;
ylabel(c, '\sigma_1 (kg m^{-3})', 'FontSize', 10); 
ylabel('Depth (m)')
%caxis([1032 1032.1]-1000)
datetick
set(gca, 'xticklabel', []);
xlim([timeSensor1(1) timeSensor1(end)])

adjust_space
pos1 = get(gca, 'pos');
adjust_space
pos2 = get(gca, 'pos');
pos2(4) = 2*pos2(4)+dy;
set(gca, 'pos', pos2)
drawnow

cbPos = get(c, 'pos');
cbPos(1) = cbPos(1)-.01;
cbPos(3) = cbPos(3)/2;
cbPos(2) = cbPos(2)+.05;
cbPos(4) = cbPos(4)-2*.05;
set(c, 'pos', cbPos)
drawnow

s2 = subplot(4,1,3);
semilogy(timeSensor1, Jb1_Vec) 
hold on
semilogy(timeSensor1, epsVec1, 'r')
datetick
set(gca, 'xticklabel', []);
xlim([timeSensor1(1) timeSensor1(end)])
%xlabel(sprintf(datestr(timeSensor(I(1)), 1)));
ylabel('J_b,\epsilon (m^2 s^{-3})')
set(gca, 'ytick', [1e-9 1e-8 1e-7 1e-6])
ylim([1e-9 2e-6])
adjust_space
% $$$ text(datenum(2012, 10, 8, 12, 27, 0), 1e-6, sprintf('J_b (2D) = %3.2g',Jb1_2D),'verticalAlignment','top')
% $$$ text(datenum(2012, 10, 8, 12, 27, 0), 5e-7, sprintf('J_b (1D) = %3.2g',Jb1_trad),'verticalAlignment','top')
% $$$ text(datenum(2012, 10, 8, 13, 10, 0), 1e-6, sprintf('eps = %3.2g',eps1_2D),'verticalAlignment','top')
% $$$ text(datenum(2012, 10, 8, 13, 10, 0), 5e-7, sprintf('eps (1D) = %3.2g',nanmean(epsVec1)),'verticalAlignment','top')
drawnow

s3 = subplot(4,1,4);
plot(timeSensor1,Jb1_Vec./epsVec1)
datetick
hold on
plot(timeSensor1,Jb1_Vec*0+.2, '--r')
ylim([0 1.5])
xlim([timeSensor1(1) timeSensor1(end)])
%xlabel(sprintf(datestr(timeSensor(I(1)), 1)));
ylabel('\Gamma');
text(timeSensor1(end)-10/1440, 1.25, sprintf('Gam = %1.2f',nanmean(Jb1_Vec./epsVec1)))
adjust_space

outfile = sprintf('Reorder2D_Jb_20121009_%02dh25.png',myTimes(ii));
print('-dpng', '-r300', outfile)
% $$$ end % end of loop on window

% $$$ figure(2)
% $$$ clf
% $$$ imagesc(timeSensor2, zVecReg, Ritp2)
% $$$ c = colorbar;
% $$$ ylabel(c, '\rho (kg m^{-3})', 'FontSize', 10); 
% $$$ xlabel(sprintf(datestr(timeSensor(I(1)), 1)));
% $$$ ylabel('Depth (m)')
% $$$ datetick

% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 4; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.04; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.14; % very right of figure
tops = 0.03; % top of figure
bots = 0.07; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %
figure(3)
clf
set(gcf, 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 16 18])
s1 = subplot(4,1,[1 2]);
imagesc(faket, fakeZ, fakeRho-1000)
c = colorbar;
ylabel(c, '\sigma_1 (kg m^{-3})', 'FontSize', 10); 
caxis([1032 1032.1])
ylabel('Depth (m)')
caxis([1032 1032.1]-1000)
datetick
set(gca, 'xticklabel', []);
xlim([faket(1) faket(end)])
adjust_space
pos1 = get(gca, 'pos');
adjust_space
pos2 = get(gca, 'pos');
pos2(4) = 2*pos2(4)+dy;
set(gca, 'pos', pos2)
drawnow

cbPos = get(c, 'pos');
cbPos(1) = cbPos(1)-.01;
cbPos(3) = cbPos(3)/2;
cbPos(2) = cbPos(2)+.05;
cbPos(4) = cbPos(4)-2*.05;
set(c, 'pos', cbPos)
drawnow

s2 = subplot(4,1,3);
semilogy(faket, Jb2_Vec) 
hold on
semilogy(faket, epsVec2, 'r')
datetick
xlim([faket(1) faket(end)])
xlabel(sprintf(datestr(timeSensor(I(1)), 1)));
ylabel('J_b,\epsilon (m^2 s^{-3})')
ylim([1e-9 2e-6])
set(gca, 'ytick', [1e-9 1e-8 1e-7 1e-6])
adjust_space
text(datenum(2012, 10, 8, 12, 27, 0), 1e-6, sprintf('J_b (2D) = %3.2g',Jb2_2D),'verticalAlignment','top')
text(datenum(2012, 10, 8, 12, 27, 0), 5e-7, sprintf('J_b (1D) = %3.2g',Jb2_trad),'verticalAlignment','top')
text(datenum(2012, 10, 8, 13, 10, 0), 1e-6, sprintf('eps (2D) = %3.2g',eps2_2D),'verticalAlignment','top')
text(datenum(2012, 10, 8, 13, 10, 0), 5e-7, sprintf('eps (1D) = %3.2g',nanmean(epsVec2)),'verticalAlignment','top')

s3 = subplot(4,1,4);
plot(faket,Jb2_Vec./epsVec2)
datetick
hold on
plot(faket,Jb2_Vec*0+.2, '--r')
ylim([0 1.5])
xlim([timeSensor1(1) timeSensor1(end)])
xlabel(sprintf(datestr(timeSensor(I(1)), 1)));
ylabel('\Gamma');
adjust_space

drawnow
print('-dpng', '-r300', 'Reorder2D_Jb_Rankine.png')


% $$$ figure(4)
% $$$ clf
% $$$ imagesc(timeSensor1, zVecReg, windowSort)
% $$$ c = colorbar;
% $$$ ylabel(c, '\rho (kg m^{-3})', 'FontSize', 10); 
% $$$ xlabel(sprintf(datestr(timeSensor(I(1)), 1)));
% $$$ caxis([1032 1032.1])
% $$$ ylabel('Depth (m)')
% $$$ caxis([1032 1032.1])
% $$$ datetick
% $$$ 
% $$$ figure(5)
% $$$ clf
% $$$ subplot(3,1,1)
% $$$ imagesc(faket, fakeZ, windowSort)
% $$$ c = colorbar;
% $$$ ylabel(c, '\rho (kg m^{-3})', 'FontSize', 10); 
% $$$ xlabel(sprintf(datestr(timeSensor(I(1)), 1)));
% $$$ caxis([1032 1032.1])
% $$$ ylabel('Depth (m)')
% $$$ caxis([1032 1032.1])
% $$$ datetick


[Jb1_2D Jb1_trad];
[Jb2_2D Jb2_trad];



keyboard


