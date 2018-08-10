
%% ----- Model ----- %%
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

uMatBc = nan(size(uMat));
vMatBc = nan(size(vMat));
sigMatBc = nan(size(sigMat));

for i = 1:length(timeVec)
    uMatBc(:,i) = uMat(:,i) - nanmean(uMat(:,i));
    vMatBc(:,i) = vMat(:,i) - nanmean(vMat(:,i));
    sigMatBc(:,i) = sigMat(:,i) - nanmean(sigMat(:,i));
end
% Reduce timeseries
I = find(timeVec >= XLIM(1) & timeVec <= XLIM(2));
sigMatBc = sigMatBc(:,I);
vMatBc = vMatBc(:,I);
uMatBc = uMatBc(:,I);
timeVec = timeVec(I);
    
figure(4)
set(gcf, 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 16 18])
clf
subplot(3,1,1)
imagesc(timeVec, zVec, sigMatBc)
datetick('x',15)
xlim(XLIM)
set(gca, 'xticklabel', []);
c = colorbar;
caxis([-.05 .05])
adjust_space
pause(.5)
c_pos = get(c, 'pos');
c_pos(3) = c_pos(3)/2;
c_pos(1) = c_pos(1)-0.02;
c_pos(2) = c_pos(2)+.02;
c_pos(4) = c_pos(4)-.04;
set(c, 'pos', c_pos);
text(XLIM(1)+1/100, 900, '\sigma_1(kg m^{-3})', 'fontWeight', 'bold')

subplot(3,1,2)
imagesc(timeVec, zuVec, vMatBc)
datetick('x',15)
xlim(XLIM)
set(gca, 'xticklabel', []);
ylabel('Depth(m)')
caxis([-.05 .05])
c = colorbar;
adjust_space
pause(.5)
c_pos = get(c, 'pos');
c_pos(3) = c_pos(3)/2;
c_pos(1) = c_pos(1)-0.02;
c_pos(2) = c_pos(2)+.02;
c_pos(4) = c_pos(4)-.04;
set(c, 'pos', c_pos);
text(XLIM(1)+1/100, 900, 'Along-isobaths (m s^{-1})', 'fontWeight', 'bold')

subplot(3,1,3)
imagesc(timeVec, zuVec, uMatBc)
datetick('x',15)
xlim(XLIM)
xlabel('One diurnal period')
caxis([-.05 .05])
c = colorbar;
adjust_space
pause(.5)
c_pos = get(c, 'pos');
c_pos(3) = c_pos(3)/2;
c_pos(1) = c_pos(1)-0.02;
c_pos(2) = c_pos(2)+.02;
c_pos(4) = c_pos(4)-.04;
text(XLIM(1)+1/100, 900, 'Cross-isobaths (m s^{-1})', 'fontWeight', 'bold')

set(c, 'pos', c_pos);

print('-dpng', '-r300', 'ModelBc.png')




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
U = -U;  % <---------------------------------- TO MATCH MODEL!!!!!!!!!!!!!!!!!!!!
V = -V;  % <---------------------------------- TO MATCH MODEL!!!!!!!!!!!!!!!!!!!!


%% Filter vel 24-h
dt = diff(abs(timeSensor(1:2))); %day
fs = 1/dt; % cpd
freq_low = 1.8; %cpd
Wn_low = freq_low/(fs/2);
[b,a] = butter(4, Wn_low);
Ufilt = nan(size(U));
Vfilt = nan(size(V));
for i = 1:size(U, 1)
    Tfilt(i,:) = filtfilt(b, a, TMat(i,:));
    Vfilt(i,:) = filtfilt(b, a, V(i,:));
    Ufilt(i,:) = filtfilt(b, a, U(i,:));    
end


%% compute bariclinic vel
uMatBc = nan(size(U));
vMatBc = nan(size(V));
TMatBc = nan(size(TMat));

for i = 1:length(timeSensor)
    TMatBc(:,i) = TMat(:,i) - nanmean(TMat(:,i));
    uMatBc(:,i) = U(:,i) - nanmean(U(:,i));
    vMatBc(:,i) = V(:,i) - nanmean(V(:,i));
end

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
   tempVec = nanmean(TMatBc(:,I),2);
   UAve(:,i) = nanmean(uMatBc(:,I),2);
   VAve(:,i) = nanmean(vMatBc(:,I),2);

   % interpolate missing sensors (in T)
   I = ~isnan(tempVec);
   tempVec = interp1(zVecReg(I), tempVec(I), zVecReg);   
   tempAve(:,i) = tempVec;
end



figure(5)
clf
set(gcf, 'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 16 18])
subplot(3,1,1)
imagesc(timeAve, zVecReg, tempAve)
datetick('x',15)
xlim([timeAve(1) timeAve(end)])
set(gca, 'xticklabel', []);
ylabel('Depth(m)')
caxis([-.05 .05])
c = colorbar;
%caxis([6.5 9])
adjust_space
pause(.5)
c_pos = get(c, 'pos');
c_pos(3) = c_pos(3)/2;
c_pos(1) = c_pos(1)-0.02;
c_pos(2) = c_pos(2)+.02;
c_pos(4) = c_pos(4)-.04;
set(c, 'pos', c_pos);
XLIM = get(gca, 'xlim');
text(XLIM(1)+1/100, 900, '\sigma_1(kg m^{-3})', 'fontWeight', 'bold')

subplot(3,1,2)
imagesc(timeAve, zVecReg, UAve)
datetick('x',15)
xlim([timeAve(1) timeAve(end)])
set(gca, 'xticklabel', []);
ylabel('Depth(m)')
caxis([-.05 .05])
c = colorbar;
adjust_space
pause(.5)
c_pos = get(c, 'pos');
c_pos(3) = c_pos(3)/2;
c_pos(1) = c_pos(1)-0.02;
c_pos(2) = c_pos(2)+.02;
c_pos(4) = c_pos(4)-.04;
set(c, 'pos', c_pos);
text(XLIM(1)+1/100, 900, 'Along-isobaths (m s^{-1})', 'fontWeight', 'bold')

subplot(3,1,3)
imagesc(timeAve, zVecReg, VAve)
xlim([timeAve(1) timeAve(end)])
datetick('x',15)
xlim([timeAve(1) timeAve(end)])
xlabel('One diurnal period')
caxis([-.05 .05])
c = colorbar;
adjust_space
pause(.5)
c_pos = get(c, 'pos');
c_pos(3) = c_pos(3)/2;
c_pos(1) = c_pos(1)-0.02;
c_pos(2) = c_pos(2)+.02;
c_pos(4) = c_pos(4)-.04;
set(c, 'pos', c_pos);
text(XLIM(1)+1/100, 900, 'Cross-isobaths (m s^{-1})', 'fontWeight', 'bold')

print('-dpng', '-r300', 'MooringBc.png')