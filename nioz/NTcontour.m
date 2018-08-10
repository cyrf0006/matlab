clear all
close all

% $$$ %% Few parameters
g = 9.81;
dz = .6;
totalDepth = 919;


%% By-pass the preceeding and load existing file
load apef_thorpe_wholeTM.mat

%% Reduce interval
I1 = find(timeSensor >= datenum(2012,10,10,20,25,0) & timeSensor <= datenum(2012,10,11,2,25,0));%Bore2
                                                                                                %I1 = find(timeSensor >= datenum(2012,10,10,14,25,0) & timeSensor <= datenum(2012,10,11,8,25,0));%Bore2
Titp = Titp(:,I1);
Ritp = Ritp1(:,I1);

timeSensor1 = timeSensor1(I1);

%% Read ADCP data
%load('/media/Seagate1TB/NIOZ/thermistordata/ROC12/roc12.mat');
load('/home/cyrf0006/research/NIOZ/RockallBank/roc12.mat');
N = [SerNmmpersec/1000]'; % in m/s now
E = [SerEmmpersec/1000]';
W = [SerVmmpersec/1000]';

zAdcp = AnDepthmm/1000; % instrument depth
timeAdcp = [datenum(SerYear+2000, SerMon, SerDay, SerHour, SerMin, SerSec)]';
%timeAdcp = timeAdcp+1/24;
% isolate good timeserie (remove if depth <10)
I = abs(zAdcp - nanmean(zAdcp))>10;
zAdcp(I) = [];
timeAdcp(I) = [];
N(:,I) = [];
E(:,I) = [];
W(:,I) = [];

% eliminate data below bottom
I = find(abs(E)>10);
E(I) = NaN;
N(I) = NaN;
W(I) = NaN;
I = sum(isnan(E),2);
II = find(I>mean(I));
Iremove = II(1):length(I);
E(Iremove,:) = [];
N(Iremove,:) = [];
W(Iremove,:) = [];

% roation may be not necessary
[U,V] = rotate_vecd(E,N,-30);

for i = 1:length(zAdcp)    
    zMatrix(:,i) = [SerBins*RDIBinSize+zAdcp(i)+RDIBin1Mid]';
end
zMatrix(Iremove, :) = [];
clear Ser* An*

I = find(timeAdcp >= timeSensor1(1) & timeAdcp <= timeSensor1(end));
U = U(:,I);
V = V(:,I);
timeAdcp = timeAdcp(I);
zMatrix = zMatrix(:,I);

% Average velocities
dzAdcp = 20;
zVecAdcp = 830:dzAdcp:950;
Uadcp = nan(length(zVecAdcp),length(timeAdcp));
Vadcp = nan(length(zVecAdcp),length(timeAdcp));

for j = 1:length(timeAdcp)
    for i = 1:length(zVecAdcp)
        I = find(zMatrix(:,j)>=zVecAdcp(i)-dzAdcp/2 & zMatrix(:,j)<=zVecAdcp(i)+dzAdcp/2);
        Uadcp(i,j) = nanmean(U(I,j));
        Vadcp(i,j) = nanmean(V(I,j));
        Wadcp(i,j) = nanmean(W(I,j));
    end
end


%% Reorder vertical profiles
disp('Reordering profiles')
windowRho = Ritp;
windowTime = timeSensor1;
H = size(windowRho,1);

% N2 calculation
[windowRhoSort, I] = sort(windowRho,1);
[dRx, dRz] = gradient(windowRhoSort);
N2bg = g./windowRhoSort.*dRz./dz;

% Compute local N2
dMat = nan(length(zVecReg), length(windowTime));
N2Mat = nan(length(zVecReg), length(windowTime));
for j = 1:length(windowTime)
    rhoVec_raw =  windowRho(:,j);
    [rhoVec_sort, I] = sort(rhoVec_raw);
    d = zVecReg - zVecReg(I);
    rho0 = nanmean(windowRho(:,j));

    dMat(:,j) = d';
    
    % local stratification
    N2local = g./rho0.*gradient(rhoVec_sort)./gradient(zVecReg');
    N2Mat(:,j) = N2local;
end
disp('  -> done!')


%% Few plots
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 2; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.1; % vert. space between subplots
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
set(gcf, 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 16 14])
V = 5:.1:10;
s1 = subplot(2,1,1);
imagesc(timeSensor1, totalDepth-zVecReg, Titp)
hold on
contour(timeSensor1, totalDepth-zVecReg, Titp, V, 'color', 'k')


%% plot ADCP

[Y, Z] = meshgrid(timeAdcp, zVecAdcp);
Z = totalDepth-Z; % swap to good z
[jmax kmax] = size(Y);
vecColor= [1 1 1];
unit = false;
scale = .3;
Vmean = nanmean(Vadcp(:));
Wmean = nanmean(Wadcp(:));
for j = 1:jmax
    for k = 1:kmax
        arrow7_no1000(Y(j,k),Z(j,k),-Vadcp(j,k),Wadcp(j,k)*0,scale,vecColor,'false');
    end
end
hold off

c = colorbar;
%ylabel(c, '\sigma_1 (kg m^{-3})', 'FontSize', 10); 
ylabel(c, 'T (^{\circ}C)', 'FontSize',10); 
ylabel('hab (m)')
%caxis([1032 1032.1]-1000)
datetick

XTIKS = get(gca, 'xtick');
time2 = time2maxV('/home/cyrf0006/research/NIOZ/RockallBank/roc12.mat',XTIKS);
time2 = round(time2*100)/100;
set(gca, 'xticklabel', time2);
set(gca, 'ydir', 'normal');
xlim([timeSensor1(1) timeSensor1(end)])
ylim([0 max(totalDepth-zVecReg)])

adjust_space
drawnow

cbPos = get(c, 'pos');
cbPos(1) = cbPos(1)-.01;
cbPos(3) = cbPos(3)/2;
cbPos(2) = cbPos(2)+.05;
cbPos(4) = cbPos(4)-2*.05;
set(c, 'pos', cbPos)
drawnow
xlabel('\phi (day)')
text(timeSensor1(end)-10/1440 , 10, 'a', 'fontSize', 12, 'fontWeight', ...
     'bold', 'color', [1 1 1])

s2 = subplot(2,1,2);
imagesc(timeSensor1, totalDepth-zVecReg, log10(N2Mat))
hold on
contour(timeSensor1, totalDepth-zVecReg, flipud(sort(Titp,1)), V, 'color', 'k')
contour(timeSensor1, totalDepth-zVecReg, Titp, V, 'color', 'k')
hold off
c = colorbar;
%ylabel(c, '\sigma_1 (kg m^{-3})', 'FontSize', 10); 
ylabel(c, 'log(N^2 / (s^{-2}))', 'FontSize', 10); 
%ylabel('Depth (m)')
%caxis([1032 1032.1]-1000)
datetick
xlim([timeSensor1(1) timeSensor1(end)])
ylim([0 max(totalDepth-zVecReg)])
set(gca, 'ydir', 'normal');

adjust_space
drawnow

cbPos = get(c, 'pos');
cbPos(1) = cbPos(1)-.01;
cbPos(3) = cbPos(3)/2;
cbPos(2) = cbPos(2)+.05;
cbPos(4) = cbPos(4)-2*.05;
set(c, 'pos', cbPos)
drawnow
caxis([-8 -3])
xlabel('10/11 October 2012')
text(timeSensor1(end)-10/1440 , 10, 'b', 'fontSize', 12, 'fontWeight', 'bold')



%outfile = sprintf('Jb-eps_%02d.eps',ii)
outfile = 'NT_toberename.eps';

print('-depsc2', outfile)

% $$$ 
% $$$ figure(2)
% $$$ clf
% $$$ imagesc(timeSensor1, totalDepth-zVecReg, dMat)
% $$$ hold on
% $$$ contour(timeSensor1, totalDepth-zVecReg, Titp, 20, 'color', 'k')
% $$$ hold off
% $$$ c = colorbar;
