%%%%%%%% ------------ Model ---------- %%%%%%%%

% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 3; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.03; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.1; % very right of figure
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
subplot(3,1,1)
contourf(timeVec, zVec, sigMat, 100, 'lineStyle', 'none')
hold on
contour(timeVec, zVec, sigMat, 50, 'color', 'k')
datetick('x',15)
xlim(XLIM)
set(gca, 'xticklabel', []);
c = colorbar;
caxis([31.9 32.1])
adjust_space
pause(.5)
c_pos = get(c, 'pos');
c_pos(3) = c_pos(3)/2;
c_pos(1) = c_pos(1)-0.02;
c_pos(2) = c_pos(2)+.02;
c_pos(4) = c_pos(4)-.04;
set(c, 'pos', c_pos);
set(gca, 'ydir', 'reverse')
text(XLIM(1)+1/100, 900, '\sigma_1(kg m^{-3})', 'fontWeight', 'bold')

subplot(3,1,2)
contourf(timeVec, zuVec, vMat, 100, 'lineStyle', 'none')
datetick('x',15)
xlim(XLIM)
set(gca, 'xticklabel', []);
ylabel('Depth(m)')
caxis([-.15 .15])
c = colorbar;
adjust_space
pause(.5)
c_pos = get(c, 'pos');
c_pos(3) = c_pos(3)/2;
c_pos(1) = c_pos(1)-0.02;
c_pos(2) = c_pos(2)+.02;
c_pos(4) = c_pos(4)-.04;
set(c, 'pos', c_pos);
set(gca, 'ydir', 'reverse')
text(XLIM(1)+1/100, 900, 'Along-isobaths (m s^{-1})', 'fontWeight', 'bold')

subplot(3,1,3)
contourf(timeVec, zuVec, uMat, 100, 'lineStyle', 'none')
datetick('x',15)
xlim(XLIM)
xlabel('One diurnal period')
caxis([-.15 .15])
c = colorbar;
set(gca, 'ydir', 'reverse')
adjust_space
pause(.5)
c_pos = get(c, 'pos');
c_pos(3) = c_pos(3)/2;
c_pos(1) = c_pos(1)-0.02;
c_pos(2) = c_pos(2)+.02;
c_pos(4) = c_pos(4)-.04;

text(XLIM(1)+1/100, 900, 'Cross-isobaths (m s^{-1})', 'fontWeight', 'bold')

set(c, 'pos', c_pos);

print('-dpng', '-r300', 'ModelBt.png')




%% $$$$$$$$$$$$$$$$----- Mooring ----- $$$$$$$$$$$$$$ %%
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 3; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.03; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.1; % very right of figure
tops = 0.03; % top of figure
bots = 0.07; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

%% Temp
wholeTS = load('~/research/NIOZ/RockallBank/apef_thorpe_wholeTM.mat');
TMat = wholeTS.Titp;
timeSensor = wholeTS.timeSensor;
zVecReg = wholeTS.zVecReg;
clear wholeTS

% to density
TSrelFile = '/home/cyrf0006/research/NIOZ/RockallBank/temp_rho_relFile_3rd.mat';
load(TSrelFile);
TMat = polyval(p, TMat)+1000;

%% ADCP
adcpFile = '/home/cyrf0006/research/NIOZ/RockallBank/roc12.mat';
[UU, E, N, W] =  velocity2hst(adcpFile, zVecReg, timeSensor);
[U, V] = rotate_vecd(E,N,-30); 
U = -U;  % <---------------------------------- TO MATCH MODEL!!!!!!!!!!!!!!!!!!!!1
V = -V;  % <---------------------------------- TO MATCH MODEL!!!!!!!!!!!!!!!!!!!!1

%% Tidal ave for mooring
% Average temperature field relative to max vel
time2 = time2maxV('/home/cyrf0006/research/NIOZ/RockallBank/roc12.mat',timeSensor); 
dtAve = .05;
timeAve = [-.5:dtAve:.5];
tempAve = nan(length(zVecReg), length(timeAve));
UAve = nan(length(zVecReg), length(timeAve));
VAve = nan(length(zVecReg), length(timeAve));
for i = 1:length(timeAve)
   I = find(time2>=timeAve(i)-dtAve & time2<=timeAve(i)+dtAve);
   tempVec = nanmean(TMat(:,I),2);
   UAve(:,i) = nanmean(U(:,I),2);
   VAve(:,i) = nanmean(V(:,I),2);

   % interpolate missing sensors (in T)
   I = ~isnan(tempVec);
   tempVec = interp1(zVecReg(I), tempVec(I), zVecReg);   
   tempAve(:,i) = tempVec;
end

figure(2)
clf
set(gcf, 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 16 18])
subplot(3,1,1)
imagesc(timeAve, zVecReg, tempAve-1000)
datetick('x',15)
xlim([timeAve(1) timeAve(end)])
set(gca, 'xticklabel', []);
caxis([31.9 32.1])
c = colorbar;
adjust_space
pause(.5)
c_pos = get(c, 'pos');
c_pos(3) = c_pos(3)/2;
c_pos(1) = c_pos(1)-0.02;
c_pos(2) = c_pos(2)+.02;
c_pos(4) = c_pos(4)-.04;
set(c, 'pos', c_pos);
XLIM2 = get(gca, 'xlim');
text(XLIM2(1)+1/100, 900, '\sigma_1(kg m^{-3})', 'fontWeight', 'bold')

subplot(3,1,2)
imagesc(timeAve, zVecReg, UAve)
datetick('x',15)
xlim([timeAve(1) timeAve(end)])
set(gca, 'xticklabel', []);
ylabel('Depth(m)')
caxis([-.15 .15])
c = colorbar;
adjust_space
pause(.5)
c_pos = get(c, 'pos');
c_pos(3) = c_pos(3)/2;
c_pos(1) = c_pos(1)-0.02;
c_pos(2) = c_pos(2)+.02;
c_pos(4) = c_pos(4)-.04;
set(c, 'pos', c_pos);
text(XLIM2(1)+1/100, 900, 'Along-isobaths (m s^{-1})', 'fontWeight', 'bold')

subplot(3,1,3)
imagesc(timeAve, zVecReg, VAve)
datetick('x',15)
xlim([timeAve(1) timeAve(end)])
xlabel('One diurnal period')
caxis([-.15 .15])
c = colorbar;
adjust_space
pause(.5)
c_pos = get(c, 'pos');
c_pos(3) = c_pos(3)/2;
c_pos(1) = c_pos(1)-0.02;
c_pos(2) = c_pos(2)+.02;
c_pos(4) = c_pos(4)-.04;
text(XLIM2(1)+1/100, 900, 'Cross-isobaths (m s^{-1})', 'fontWeight', 'bold')

set(c, 'pos', c_pos);

print('-dpng', '-r300', 'MooringBt.png')


%% Sig only (WTD pres)
% *********************** Adjust_space.m ************************ %
% Fields required by the function adjust_space.m. Please fill every
% of the following and call "adjust_space" in the script whenever
% you want. Do not touch four last fields
ncol = 1; % no. subplot column
nrow = 1; % no. subplot row
dx = 0.03 ; % horiz. space between subplots
dy = 0.03; % vert. space between subplots
lefs = 0.1; % very left of figure
rigs = 0.12; % very right of figure
tops = 0.03; % top of figure
bots = 0.1; % bottom of figure
figw = (1-(lefs+rigs+(ncol-1)*dx))/ncol;
figh = (1-(tops+bots+(nrow-1)*dy))/nrow;
count_col = 1;
count_row = 1;
% *************************************************************** %

figure(3)
clf
set(gcf, 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 14 6])
contourf(timeVec, zVec, sigMat, 100, 'lineStyle', 'none')
hold on
contour(timeVec, zVec, sigMat, 50, 'color', 'k')
datetick('x',15)
xlim(XLIM)
c = colorbar;
caxis([31.9 32.1])
adjust_space
pause(.5)
c_pos = get(c, 'pos');
c_pos(3) = c_pos(3)/2;
c_pos(1) = c_pos(1)-0.02;
c_pos(2) = c_pos(2)+.02;
c_pos(4) = c_pos(4)-.04;
set(c, 'pos', c_pos);
set(gca, 'ydir', 'reverse')
text(XLIM(1)+1/100, 910, '\sigma_1(kg m^{-3})', 'fontWeight', 'bold')


xlabel('One diurnal period')
print('-dpng', '-r300', 'ModelBore.png')
